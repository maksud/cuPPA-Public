/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "cuPPA.cuh"
#include <helper_cuda.h>

#if USE_DEFAULT_CURAND
__device__ bool re_execute_copy_model_rrp(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandState *state)
#else
__device__ bool re_execute_copy_model_rrp(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandStateMRG32k3a *state)
#endif
{
#if DEBUG_cuPPA_RRP
    int nRetry = 0;
    while (nRetry < 100)
#else
    while (true)
#endif
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
#if DEBUG_cuPPA_RRP
        nRetry++;
#endif
    }
#if DEBUG_cuPPA_RRP
    cuDBG(2, "Which Node: %d\n", u);
    V->printInfo();
    cuDBG(2, "MAX RETRY FAILED\n");
#endif
    return false;
}

__global__ void initial_graph_rrp(EdgeIndexType M, NodeType *G)
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    cuDBG0(2, "nThreads: %d, %d\n", nThreads, M);

    //m disjoint graphs
    for (int i = id; i < M; i += nThreads)
    {
        for (int j = 0; j < M; j++)
        {
            G[i * M + j] = i;
        }
    }

    //m-th node
    if ((M % nThreads) == id) //Select a Specific Thread, also can be arbitrary!
    {
        for (int j = 0; j < M; j++)
        {
            G[M * M + j] = j;
        }
    }

    cuDBG0(2, "Thread %d: Initial Graph Generated!\n", id);
}

__global__ void initialize_vertex_structures_rrp(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F)
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    cuDBG0(2, "nThreads: %d, %d\n", nThreads, M);

    for (int i = id; i < B; i += nThreads)
    {
        F[i]._M = M;
        F[i]._count = 0;
        F[i]._root = 0;

        F[i]._A = &G[(u + i) * M];
        F[i]._leftChild = &children[(2 * i) * M];
        F[i]._rightChild = &children[(2 * i + 1) * M];

        for (int j = 0; j < M; j++)
        {
            F[i]._A[j] = -1;
            F[i]._leftChild[j] = NULL_PTR;
            F[i]._rightChild[j] = NULL_PTR;
        }
    }
    cuDBG0(2, "Thread %d: Vertex Setup Completed!\n", id);
}

__global__ void initialize_vertex_structures_rrp_fast(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F)
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    cuDBG0(2, "nThreads: %d, %d\n", nThreads, M);

    for (int i = id; i < B; i += nThreads)
    {
        F[i]._M = M;
        F[i]._count = 0;
        F[i]._root = 0;

        F[i]._A = &G[(u + i) * M];

//        for (int j = 0; j < M; j++)
//        {
//            F[i]._A[j] = -1;
//            F[i]._leftChild[j] = NULL_PTR;
//            F[i]._rightChild[j] = NULL_PTR;
//        }
    }
    cuDBG0(2, "Thread %d: Vertex Setup Completed!\n", id);
}

__global__ void initialize_vertex_structures_rrp_2(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    cuDBG0(2, "nThreads: %d, %d\n", blockDim.x * gridDim.x, M);

    if (id < B)
    {
        F[id]._M = M;
        F[id]._count = 0;
        F[id]._root = 0;

        F[id]._A = &G[(u + id) * M];
        F[id]._leftChild = &children[(2 * id) * M];
        F[id]._rightChild = &children[(2 * id + 1) * M];

        for (int j = 0; j < M; j++)
        {
            F[id]._A[j] = -1;
            F[id]._leftChild[j] = NULL_PTR;
            F[id]._rightChild[j] = NULL_PTR;
        }
    }
    cuDBG0(2, "Thread %d: Vertex Setup Completed!\n", id);
}

__global__ void validate1_graph_rrp_by_batch(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F)
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = id; i < B; i += nThreads)
    {
        for (int j = 0; j < M; j++)
        {
            if (G[(u + i) * M + j] != -1)
            {
                cuDBG(2, "Validate1 %ld,%d\n", u + i, j);
            }
        }
    }
}

__global__ void validate2_graph_rrp_by_batch(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F)
{
    int nThreads = blockDim.x * gridDim.x;
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = id; i < B; i += nThreads)
    {
        for (int j = 0; j < M; j++)
        {
            if (G[(u + i) * M + j] < 0)
            {
                cuDBG(2, "############# ############# ############# ############# ############# Validate2 %lld,%d\t%lld\n", u + i, j, G[(u + i) * M + j]);
            }
        }
    }
}
