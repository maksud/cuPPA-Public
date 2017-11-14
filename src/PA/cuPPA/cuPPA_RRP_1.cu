/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "cuPPA.cuh"
#include <helper_cuda.h>

__global__ void find_max_queue_size(int X[])
{
    __shared__ int WQ[6 * 1024];

    WQ[threadIdx.x] = X[threadIdx.x];
    __syncthreads();
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (threadIdx.x < s)
        {
            WQ[threadIdx.x] = max(WQ[threadIdx.x], WQ[threadIdx.x + s]);
        }
        __syncthreads();
    }

    if (threadIdx.x == 0)
        X[0] = WQ[0];
}

#if USE_DEFAULT_CURAND
__global__ void generate_graph_rrp_by_batch(const NodeType startingV, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandState *state, int *MaxQSize)
#else
__global__ void generate_graph_rrp_by_batch(const NodeType startingV, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state, int *MaxQSize)
#endif
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    cuDBG0(2, "nThreads: %d, %ld\t%d\t%f\n", nThreads, B, M, p);
    cuDBG0(2, "Processing: Vertex %ld to %ld\n", startingV, startingV + B - 1);

    const NodeType lb = startingV + id;
    const NodeType ub = startingV + B;

#if USE_DEFAULT_CURAND
    curandState localState = state[id]; //Save Random Number State to the global memory
#else
    curandStateMRG32k3a localState = state[id]; //Save Random Number State to the global memory
#endif

    const int MAX_QUEUE_SIZE = (48 * 1024 / sizeof(EdgeType)); // 48KB Shared Memory can store MAX_QUEUE_SIZE items

    __shared__ EdgeType WQ[MAX_QUEUE_SIZE]; //3K*8=24KB
    cuDBG0(2, "sizeof(WQ): %ld\n", sizeof(WQ));
    cuDBG0(2, "MAX_QUEUE_SIZE: %d\n", MAX_QUEUE_SIZE);
    cuDBG0(2, "blockDim.x: %d\n", blockDim.x);

    const int MAX_SIZE_PER_THREAD = MAX_QUEUE_SIZE / blockDim.x;
    cuDBG0(2, "MAX_SIZE_PER_THREAD: %d\n", MAX_SIZE_PER_THREAD);

    int wqIndexCurrent = 0;
    int MAX_WQ_SIZE = 0;
    EdgeType* WQCurrent;

    WQCurrent = &WQ[threadIdx.x * MAX_SIZE_PER_THREAD];

    for (NodeType u = lb; u < ub; u += nThreads)
    {

        for (EdgeIndexType i = 0; i < M;)
        {
            NodeType v = curand(&localState) % u;
            double r = curand_uniform(&localState);

            if (r < p)
            {
                //F_i(u) <- v
                if (F[u - startingV].insertItem(i, v))
                {
                    i++; //Success, next outgoing edge
                }
            }
            else
            {
                //F_i(u) <- F_k(v)

                EdgeIndexType k = curand(&localState) % M;

                NodeType F_k_v = G[v * M + k];

                if (F_k_v >= 0) //Try to resolve the items first.
                {
                    //F_k(v) is Known, Try to insert or retry
                    if (F[u - startingV].insertItem(i, F_k_v))
                    {
                        i++; //Success, next outgoing edge
                    }
                }
                else
                {
                    //F_k(v) is not Known, Insert a placeholder if possible. Otherwise retry!
                    EdgeType E_Out = u * M + i;
                    EdgeType E_In = -(v * M + k); // -ve indicates Incomplete!

                    if (wqIndexCurrent < MAX_SIZE_PER_THREAD)
                    {
                        if (F[u - startingV].insertItem(i, E_In))
                        {
                            cuDBG(2, "%d: Inserted %ld @ %d\n", id, E_Out, wqIndexCurrent);
                            WQCurrent[wqIndexCurrent++] = E_Out;
                            i++;
                        }
                    }
                    else
                    {
                        //Not enough space on queue. Re-run Copy Model
//                        re_execute_copy_model_rrp(G, &F[u - startingV], M, p, w, i, &localState); //
//                        i++;
                        cuDBG(2, "OVERFLOW\n");
                    }
                }
            }
        }
//        cuDBG(5, "%d: Worked on: %ld,%d\n", id, u, F[u - startingV].count());
    }

//    cuDBG(2, "%d: Max WQ Size: %d\n", id, wqIndexCurrent);

    MAX_WQ_SIZE = wqIndexCurrent;

    int wqIndexNew = 0;
    int nResolved;
    int nReRunCopyModel;
    int nStep = 0;

    while (wqIndexCurrent > 0)
    {
        {
            nStep++;
//            cuDBG(2, "STEP %d: ******************************************\n", nStep);
//            cuDBG(4, "%d B1: Max WQ Size: NEW: %d\tCURRENT: %d\t nStep: %d \n", id, wqIndexNew, wqIndexCurrent, nStep);
        }

        nResolved = 0;
        nReRunCopyModel = 0;

        for (int q_i = 0; q_i < wqIndexCurrent; q_i++)
        {
            EdgeType e = WQCurrent[q_i];
            NodeType w = e / M;
            EdgeIndexType i = e % M;

            EdgeType eOut = -G[w * M + i];
            NodeType w2 = eOut / M;
            EdgeIndexType i2 = eOut % M;

            //Mutually Exclusive Operations!
            if (G[w2 * M + i2] >= 0)
            {
                if (F[w - startingV].updateItem(F[w - startingV].Adj(i), G[w2 * M + i2]) >= 0)
                {
                    nResolved++;
                }
                else
                {
                    //Can not resolve as duplicate edge will form. Re-run Copy Model
                    re_execute_copy_model_rrp(G, &F[w - startingV], M, p, w, i, &localState);
                    nReRunCopyModel++;
                }
            }
            else
            {
                WQCurrent[wqIndexNew++] = e;

            }
        }

//        cuDBG(2, "%d B2: Max WQ Size: NEW: %d\t CURRENT: %d\tRESOLVED: %d\tRERUN: %d\n", id, wqIndexNew, wqIndexCurrent, nResolved, nReRunCopyModel);

        {
            wqIndexCurrent = wqIndexNew;
            wqIndexNew = 0;

            nReRunCopyModel = 0;
            nResolved = 0;
        }
    }
    __syncthreads();

//    cuDBG(2, "Step: %d\n", nStep);
    state[id] = localState; //Restore Random Number State to the global memory

    MaxQSize[blockIdx.x] = 0;
    WQ[threadIdx.x] = (int) MAX_WQ_SIZE;
    __syncthreads();
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (threadIdx.x < s)
        {
            WQ[threadIdx.x] = max((int) WQ[threadIdx.x], (int) WQ[threadIdx.x + s]);
        }
        __syncthreads();
    }
    __syncthreads();

    if (threadIdx.x == 0)
    {
        MaxQSize[blockIdx.x] = (int) WQ[0];
//        printf("%d: MaxQSize = %d\n", blockIdx.x, (int) WQ[0]);
    }
}
