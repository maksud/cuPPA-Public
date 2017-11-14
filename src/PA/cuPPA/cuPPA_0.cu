/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "cuPPA.cuh"
#include <helper_cuda.h>

#if USE_DEFAULT_CURAND
__device__ bool re_execute_copy_model(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandState *state)
#else
__device__ bool re_execute_copy_model(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandStateMRG32k3a *state)
#endif
{
    int nRetry = 0;
//    while (true)
    while (nRetry < 100)
    {
        NodeType v = curand(state) % u;
        double r = curand_uniform(state);
        if (r < p)
        {
            if (V->updateItem(V->Adj(i), v) >= 0)
                return true;
        }
        else
        {
            EdgeIndexType k = curand(state) % M;
            NodeType F_k_v = G[v * M + k];
            if (F_k_v >= 0) //Try to resolve the items first.
            {
                if (V->updateItem(V->Adj(i), F_k_v) >= 0)
                    return true;
            }
        }
        nRetry++;
    }

    printf("Which Node: %d\n", u);
    V->printInfo();
    printf("MAX RETRY FAILED\n");
    return false;
}

__global__ void initialize_vertex_structures(NodeType nVertices, EdgeIndexType M, NodeType *A, EdgeIndexType *children, cuBSTVertex *vertices)
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id == 0)
        cuDBG(2, "nThreads: %d, %d\n", nThreads, M);

    for (int i = id; i < nVertices; i += nThreads)
    {
        vertices[i]._M = M;
        vertices[i]._count = 0;

//        printf("Thread %d: %d\n", id, i);

        vertices[i]._A = &A[i * M];
        vertices[i]._leftChild = &children[(2 * i) * M];
        vertices[i]._rightChild = &children[(2 * i + 1) * M];

        for (int j = 0; j < M; j++)
        {
            vertices[i]._A[j] = -1;
            vertices[i]._leftChild[j] = NULL_PTR;
            vertices[i]._rightChild[j] = NULL_PTR;
        }
    }
    if (id == 0)
        cuDBG(2, "Thread %d: Vertex Setup Completed!\n", id);
}

__global__ void initial_graph(EdgeIndexType M, cuBSTVertex *vertices)
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id == 0)
        cuDBG(2, "nThreads: %d, %d\n", nThreads, M);

    //m disjoint graphs
    for (int i = id; i < M; i += nThreads)
    {
//        printf("disconnected nodes: Thread %d: %d\n", id, i);
        for (int j = 0; j < M; j++)
        {
            vertices[i].Adj(j, i);
        }
    }

    __syncthreads();

    //m-th node
    if ((M % nThreads) == id) //Select a Specific Thread, also can be arbitrary!
    {
//        printf("m-th node: Thread %d\n", id);
        for (int j = 0; j < M; j++)
        {
            vertices[M].insertItem(j, j);
        }
    }

    if (id == 0)
        cuDBG(2, "Thread %d: Initial Graph Generated!\n", id);
}

#if USE_DEFAULT_CURAND
__global__ void generate_graph_race(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandState *state)
#else
__global__ void generate_graph_race(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state)
#endif
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    cuDBG(2, "nThreads: %d, %ld\t%d\t%f\n", nThreads, nVertices, M, p);

    NodeType lb = ((NodeType) (M / nThreads)) * nThreads + id;

    if (lb < M + 1)
        lb += nThreads;

#if USE_DEFAULT_CURAND
    curandState localState = state[id]; //Save Random Number State to the global memory
#else
    curandStateMRG32k3a localState = state[id]; //Save Random Number State to the global memory
#endif

    const int MAX_SIZE = 3000;

    __shared__ int wqIndexCurrent;
    __shared__ int wqIndexNew;

    wqIndexCurrent = 0;
    wqIndexNew = 0;

    __shared__ EdgeType WQ_1[MAX_SIZE]; //3K*8=24KB
    __shared__ EdgeType WQ_2[MAX_SIZE]; //3K*8=24KB

    __shared__ EdgeType* WQCurrent;
    __shared__ EdgeType* WQNew;
    __shared__ EdgeType* WQTmp;

    WQNew = WQ_1;
    WQCurrent = WQ_2;

    __syncthreads();

    cuDBG(2, "WQ Size: %ld\t%ld\n", sizeof(WQ_1), sizeof(WQ_2));

    for (NodeType u = lb; u < nVertices; u += nThreads)
    {
        //printf("Working Node: Thread %d: %ld\n", id, u);

        for (EdgeIndexType i = 0; i < M;)
        {
            NodeType v = curand(&localState) % u;

            double r = curand_uniform(&localState);

            if (r < p)
            {
                //F_i(u) <- v
                if (F[u].insertItem(i, v))
                    i++;
            }
            else
            {
                EdgeIndexType k = curand(&localState) % M;

                //F_i(u) <- F_k(v)

                NodeType F_k_v = G[v * M + k];

                if (F_k_v >= 0) //Try to resolve the items first.
                {
                    //F_k(v) is Known, Try to insert or retry
                    if (F[u].insertItem(i, F_k_v))
                    {
                        i++;
                    }
                }
                else
                {
                    //F_k(v) is not Known, Insert a placeholder if possible. Otherwise retry!
                    EdgeType E_Out = u * M + i;
                    EdgeType E_In = -(v * M + k); // -ve indicates Incomplete!

                    if (wqIndexCurrent < MAX_SIZE)
                    {
                        if (F[u].insertItem(i, E_In))
                        {
//                            printf("Putting: %ld into %ld\n", E_In, E_Out);

                            int storeIndex = atomicAdd(&wqIndexCurrent, 1);
                            WQCurrent[storeIndex] = E_Out;
                            i++;
                        }
                    }
                }
            }
        }
    }

    if (threadIdx.x == 0)
    {
//        F[19].printInfo();
//        printf("DEBUG:F[19,4] %ld, %d\n", F[19].Adj(4), F[19]._count);

        cuDBG(2, "Max WQ Size: %d\n", wqIndexCurrent);
    }
    __shared__ int nResolved;
    __shared__ int nReRunCopyModel;
    __shared__ int nStep;
    nStep = 0;

    while (wqIndexCurrent > 0 && nStep < 10)
    {
        if (threadIdx.x == 0)
        {
            nStep++;
            cuDBG(2, "STEP %d: ******************************************\n", nStep);

            printf("%d B1: Max WQ Size: NEW: %d\tCURRENT: %d\n", id, wqIndexNew, wqIndexCurrent);

            printf("BD\t%d\n", blockDim.x);
        }

        nResolved = 0;
        nReRunCopyModel = 0;

//        __syncthreads();
        for (int q_i = threadIdx.x; q_i < wqIndexCurrent; q_i += blockDim.x)
        {
            printf("Thread: %d\t%d\n", threadIdx.x, q_i);

            //Thread id works on items: id, id+nThreads, ..., wqIndex in a round-robin fasion.
            EdgeType e = WQCurrent[q_i];
            NodeType w = e / M;
            EdgeIndexType i = e % M;

            EdgeType eOut = -F[w].Adj(i);
            NodeType w2 = eOut / M;
            EdgeIndexType i2 = eOut % M;

            printf("Woorking on %ld (%ld,%d)\n", e, w, i);

            //Mutually Exclusive Operations!
            if (F[w2].Adj(i2) >= 0)
            {
                if (w == 229)
                {
                    printf("%d TEST 1\n", threadIdx.x);
                    F[w].printInfo();
                }

                if (F[w].updateItem(F[w].Adj(i), F[w2].Adj(i2)) >= 0)
                {
                    if (w == 229)
                    {
                        printf("Updating Debug: This Should not Happen! %ld\t%ld\n", F[w].Adj(i), F[w2].Adj(i2));

                        printf("TEST 3\n");
                        F[w].printInfo();
                    }
                    printf("Resolved, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                    atomicAdd(&nResolved, 1);
                }
                else
                {

                    if (w == 229)
                    {
                        printf("TEST 2\n");
                        F[w].printInfo();
                    }
//                    F[w].Adj(i, F[w2].Adj(i2)); //Not Checked for Duplicates YET!

                    if (re_execute_copy_model(G, &F[w], M, p, w, i, &localState))
                    {
//                        cuDBG(2, "Second-Resolved-PASS, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                        atomicAdd(&nReRunCopyModel, 1);
                    }
                    else
                    {
//                        F[w].Adj(i, F[w2].Adj(i2)); //Not Checked for Duplicates YET!
//                        cuDBG(2, "Second-Resolved-FAIL, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                    }
                }
            }
            else
            {
                //Insert into WQNew
//                printf("Un-Resolved ,%ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));

                if (wqIndexNew < MAX_SIZE)
                {
//                    printf("Re-Putting: %ld into %ld\n", eOut, e);

                    int storeIndex = atomicAdd(&wqIndexNew, 1);
                    WQNew[storeIndex] = e;
                }
            }

            __syncthreads();
        }

        __syncthreads();

        if (threadIdx.x == 0)
        {
            cuDBG(2, "%d B2: Max WQ Size: NEW: %d\t CURRENT: %d\tRESOLVED: %d\tRERUN: %d\n", id, wqIndexNew, wqIndexCurrent, nResolved, nReRunCopyModel);
        }

        if (threadIdx.x == 0)
        {
            wqIndexCurrent = wqIndexNew;
            wqIndexNew = 0;

            WQTmp = WQNew;
            WQNew = WQCurrent;
            WQCurrent = WQTmp;

            nReRunCopyModel = 0;
            nResolved = 0;
        }
        __syncthreads();
    }

    state[id] = localState; //Restore Random Number State to the global memory

    cuDBG(2, "Thread %d: Graph Generated!\n", id);

    if (id == 0)
    {
//        F[16].printInfo();
    }
}

#if USE_DEFAULT_CURAND
__global__ void generate_graph_mutex(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandState *state)
#else
__global__ void generate_graph_mutex(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state)
#endif
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    cuDBG(2, "nThreads: %d, %ld\t%d\t%f\n", nThreads, nVertices, M, p);

    NodeType lb = ((NodeType) (M / nThreads)) * nThreads + id;

    if (lb < M + 1)
        lb += nThreads;

#if USE_DEFAULT_CURAND
    curandState localState = state[id]; //Save Random Number State to the global memory
#else
    curandStateMRG32k3a localState = state[id]; //Save Random Number State to the global memory
#endif

    const int MAX_SIZE = 3000;

    __shared__ int wqIndexCurrent;
    __shared__ int wqIndexNew;

    wqIndexCurrent = 0;
    wqIndexNew = 0;

    __shared__ EdgeType WQ_1[MAX_SIZE]; //3K*8=24KB
    __shared__ EdgeType WQ_2[MAX_SIZE]; //3K*8=24KB

    __shared__ EdgeType* WQCurrent;
    __shared__ EdgeType* WQNew;
    __shared__ EdgeType* WQTmp;

    WQNew = WQ_1;
    WQCurrent = WQ_2;

    __syncthreads();

    cuDBG(2, "WQ Size: %ld\t%ld\n", sizeof(WQ_1), sizeof(WQ_2));

    for (NodeType u = lb; u < nVertices; u += nThreads)
    {
        //printf("Working Node: Thread %d: %ld\n", id, u);

        for (EdgeIndexType i = 0; i < M;)
        {
            NodeType v = curand(&localState) % u;

            double r = curand_uniform(&localState);

            if (r < p)
            {
                //F_i(u) <- v
                if (F[u].insertItem(i, v))
                    i++;
            }
            else
            {
                EdgeIndexType k = curand(&localState) % M;

                //F_i(u) <- F_k(v)

                NodeType F_k_v = G[v * M + k];

                if (F_k_v >= 0) //Try to resolve the items first.
                {
                    //F_k(v) is Known, Try to insert or retry
                    if (F[u].insertItem(i, F_k_v))
                    {
                        i++;
                    }
                }
                else
                {
                    //F_k(v) is not Known, Insert a placeholder if possible. Otherwise retry!
                    EdgeType E_Out = u * M + i;
                    EdgeType E_In = -(v * M + k); // -ve indicates Incomplete!

                    if (wqIndexCurrent < MAX_SIZE)
                    {
                        if (F[u].insertItem(i, E_In))
                        {
//                            printf("Putting: %ld into %ld\n", E_In, E_Out);

                            int storeIndex = atomicAdd(&wqIndexCurrent, 1);
                            WQCurrent[storeIndex] = E_Out;
                            i++;
                        }
                    }
                }
            }
        }
    }

    if (threadIdx.x == 0)
    {
//        F[19].printInfo();
//        printf("DEBUG:F[19,4] %ld, %d\n", F[19].Adj(4), F[19]._count);

        cuDBG(2, "Max WQ Size: %d\n", wqIndexCurrent);
    }
    __shared__ int nResolved;
    __shared__ int nReRunCopyModel;
    __shared__ int nStep;
    nStep = 0;

    if (id == 0)
        F[229].printInfo();

    while (wqIndexCurrent > 0 && nStep < 10)
    {
        if (threadIdx.x == 0)
        {
            nStep++;
            cuDBG(2, "STEP %d: ******************************************\n", nStep);

            printf("%d B1: Max WQ Size: NEW: %d\tCURRENT: %d\n", id, wqIndexNew, wqIndexCurrent);

            printf("BD\t%d\n", blockDim.x);
        }

        nResolved = 0;
        nReRunCopyModel = 0;

//        __syncthreads();
        for (int q_i = threadIdx.x; q_i < wqIndexCurrent; q_i += blockDim.x)
        {
            printf("Thread: %d\t%d\n", threadIdx.x, q_i);

            //Thread id works on items: id, id+nThreads, ..., wqIndex in a round-robin fasion.
            EdgeType e = WQCurrent[q_i];
            NodeType w = e / M;
            EdgeIndexType i = e % M;

            EdgeType eOut = -F[w].Adj(i);
            NodeType w2 = eOut / M;
            EdgeIndexType i2 = eOut % M;

            printf("Woorking on %ld (%ld,%d)\n", e, w, i);

            //Mutually Exclusive Operations!
            if (F[w2].Adj(i2) >= 0)
            {
                if (w == 229)
                {
                    printf("%d TEST 1\n", threadIdx.x);
                    F[w].printInfo();
                }

                if (F[w].updateItem(F[w].Adj(i), F[w2].Adj(i2)) >= 0)
                {
                    if (w == 229)
                    {
                        printf("Updating Debug: This Should not Happen! %ld\t%ld\n", F[w].Adj(i), F[w2].Adj(i2));

                        printf("TEST 3\n");
                        F[w].printInfo();
                    }
                    printf("Resolved, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                    atomicAdd(&nResolved, 1);
                }
                else
                {

                    if (w == 229)
                    {
                        printf("TEST 2\n");
                        F[w].printInfo();
                    }
//                    F[w].Adj(i, F[w2].Adj(i2)); //Not Checked for Duplicates YET!

                    if (re_execute_copy_model(G, &F[w], M, p, w, i, &localState))
                    {
//                        cuDBG(2, "Second-Resolved-PASS, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                        atomicAdd(&nReRunCopyModel, 1);
                    }
                    else
                    {
//                        F[w].Adj(i, F[w2].Adj(i2)); //Not Checked for Duplicates YET!
//                        cuDBG(2, "Second-Resolved-FAIL, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                    }
                }
            }
            else
            {
                //Insert into WQNew
//                printf("Un-Resolved ,%ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));

                if (wqIndexNew < MAX_SIZE)
                {
//                    printf("Re-Putting: %ld into %ld\n", eOut, e);

                    int storeIndex = atomicAdd(&wqIndexNew, 1);
                    WQNew[storeIndex] = e;
                }
            }

            __syncthreads();
        }

        __syncthreads();

        if (threadIdx.x == 0)
        {
            cuDBG(2, "%d B2: Max WQ Size: NEW: %d\t CURRENT: %d\tRESOLVED: %d\tRERUN: %d\n", id, wqIndexNew, wqIndexCurrent, nResolved, nReRunCopyModel);
        }

        if (threadIdx.x == 0)
        {
            wqIndexCurrent = wqIndexNew;
            wqIndexNew = 0;

            WQTmp = WQNew;
            WQNew = WQCurrent;
            WQCurrent = WQTmp;

            nReRunCopyModel = 0;
            nResolved = 0;
        }
        __syncthreads();
    }

    state[id] = localState; //Restore Random Number State to the global memory

    cuDBG(2, "Thread %d: Graph Generated!\n", id);
}

#if USE_DEFAULT_CURAND
__global__ void generate_graph_avoid(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandState *state)
#else
__global__ void generate_graph_avoid(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state)
#endif
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id == 0)
        cuDBG(2, "nThreads: %d, %ld\t%d\t%f\n", nThreads, nVertices, M, p);

    NodeType lb = ((NodeType) (M / nThreads)) * nThreads + id;

    if (lb < M + 1)
        lb += nThreads;

#if USE_DEFAULT_CURAND
    curandState localState = state[id]; //Save Random Number State to the global memory
#else
    curandStateMRG32k3a localState = state[id]; //Save Random Number State to the global memory
#endif

    const int MAX_QUEUE_SIZE = (48 * 1024 / sizeof(EdgeType)); // 48KB Shared Memory can store MAX_QUEUE_SIZE items
    __shared__ EdgeType WQ[MAX_QUEUE_SIZE]; //3K*8=24KB

    cuDBG(2, "MAX_QUEUE_SIZE: %d\n", MAX_QUEUE_SIZE);
    cuDBG(2, "blockDim.x: %d\n", blockDim.x);
    const int MAX_SIZE_PER_THREAD = MAX_QUEUE_SIZE / blockDim.x;
    cuDBG(2, "MAX_SIZE_PER_THREAD: %d\n", MAX_SIZE_PER_THREAD);

    int wqIndexCurrent = 0;
    int wqIndexNew = 0;

    EdgeType* WQCurrent;

    WQCurrent = &WQ[threadIdx.x * MAX_SIZE_PER_THREAD];

    if (id == 0)
        cuDBG(2, "WQ Size: %ld\n", sizeof(WQ));

    for (NodeType u = lb; u < nVertices; u += nThreads)
    {
        //printf("Working Node: Thread %d: %ld\n", id, u);
        for (EdgeIndexType i = 0; i < M;)
        {
            NodeType v = curand(&localState) % u;

            double r = curand_uniform(&localState);

            if (r < p)
            {
                //F_i(u) <- v
                if (F[u].insertItem(i, v))
                    i++;
            }
            else
            {
                EdgeIndexType k = curand(&localState) % M;

                //F_i(u) <- F_k(v)

                NodeType F_k_v = G[v * M + k];

                if (F_k_v >= 0) //Try to resolve the items first.
                {
                    //F_k(v) is Known, Try to insert or retry
                    if (F[u].insertItem(i, F_k_v))
                    {
                        i++;
                    }
                }
                else
                {
                    //F_k(v) is not Known, Insert a placeholder if possible. Otherwise retry!
                    EdgeType E_Out = u * M + i;
                    EdgeType E_In = -(v * M + k); // -ve indicates Incomplete!

                    if (wqIndexCurrent < MAX_SIZE_PER_THREAD)
                    {
                        if (F[u].insertItem(i, E_In))
                        {
//                            if (threadIdx.x == 0)
//                            {
//                                printf("%d: Putting: %ld into %ld\n", threadIdx.x, E_In, E_Out);
//                            }
                            WQCurrent[wqIndexCurrent++] = E_Out;
                            i++;
                        }
                    }
                    else
                    {
                        printf("OVERFLOW\n");
                    }
                }
            }
        }
    }

    __syncthreads();

    cuDBG(2, "%d: Max WQ Size: %d\n", threadIdx.x, wqIndexCurrent);

    int nResolved;
    int nReRunCopyModel;
    int nStep = 0;
//
//    for (int i = 0; i < wqIndexCurrent; i++)
//    {
//        if (threadIdx.x == 0)
//            printf("Woorking on %ld\n", WQCurrent[i]);
//    }
//
//    return;

    while (wqIndexCurrent > 0 && nStep < 1024)
    {
        {
            nStep++;
            cuDBG(2, "STEP %d: ******************************************\n", nStep);
            printf("%d B1: Max WQ Size: NEW: %d\tCURRENT: %d\t nStep: %d \n", id, wqIndexNew, wqIndexCurrent, nStep);
        }

        nResolved = 0;
        nReRunCopyModel = 0;

//        __syncthreads();
        for (int q_i = 0; q_i < wqIndexCurrent; q_i++)
        {
#if 0
            printf("Thread: %d\t%d\n", threadIdx.x, q_i);
#endif

            EdgeType e = WQCurrent[q_i];
            NodeType w = e / M;
            EdgeIndexType i = e % M;

            EdgeType eOut = -F[w].Adj(i);
            NodeType w2 = eOut / M;
            EdgeIndexType i2 = eOut % M;

#if 0
            if (threadIdx.x == 0)
            {
                printf("Woorking on %ld (%ld,%d)\n", e, w, i);
            }
#endif
            //Mutually Exclusive Operations!
            if (F[w2].Adj(i2) >= 0)
            {
                if (F[w].updateItem(F[w].Adj(i), F[w2].Adj(i2)) >= 0)
                {
//                    printf("Resolved, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                    nResolved++;
                }
                else
                {
                    if (re_execute_copy_model(G, &F[w], M, p, w, i, &localState))
                    {
//                        cuDBG(2, "Second-Resolved-PASS, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                        nReRunCopyModel++;
                    }
                    else //Will NEVER OCCUR!!!
                    {
                        F[w].Adj(i, F[w2].Adj(i2)); //Not Checked for Duplicates YET!
//                        cuDBG(2, "Second-Resolved-FAIL, %ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));
                    }
                }
            }
            else
            {
                //Insert into WQNew
//                printf("Un-Resolved ,%ld, (F[,%ld,%d,]) -> ,%ld, (F[,%ld,%d,]) = ,%ld\n", e, w, i, eOut, w2, i2, F[w2].Adj(i2));

                if (wqIndexNew < MAX_SIZE_PER_THREAD)
                {
//                    printf("Re-Putting: %ld into %ld\n", eOut, e);

                    WQCurrent[wqIndexNew++] = e;
                }
            }
        }

//        {
//            cuDBG(2, "%d B2: Max WQ Size: NEW: %d\t CURRENT: %d\tRESOLVED: %d\tRERUN: %d\n", id, wqIndexNew, wqIndexCurrent, nResolved, nReRunCopyModel);
//        }

        {
            wqIndexCurrent = wqIndexNew;
            wqIndexNew = 0;

            nReRunCopyModel = 0;
            nResolved = 0;
        }
        __syncthreads();
    }

    state[id] = localState; //Restore Random Number State to the global memory

    cuDBG(2, "Thread %d: Graph Generated!\n", id);
}
