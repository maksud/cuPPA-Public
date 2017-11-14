/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "cuPPA.cuh"
#include <helper_cuda.h>

__global__ void generate_graph_serial(const NodeType startingV, const NodeType B, const EdgeIndexType M, const double p, NodeType *G, curandStateMRG32k3a *state)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    cuDBG0(2, "nThreads: %d, %ld\t%d\t%f\n", blockDim.x * gridDim.x, B, M, p);
    cuDBG0(2, "Processing: Vertex %ld to %ld\n", startingV, startingV + B - 1);

    if (id == 0)
    {

        const NodeType lb = startingV + id;
        const NodeType ub = startingV + B;

        curandStateMRG32k3a localState = state[id]; //Save Random Number State to the global memory

        cuBSTVertex V(M, NULL);

        for (NodeType u = lb; u < ub; u++)
        {
            V._A = &G[lb * M];
            V.clear();

            for (EdgeIndexType i = 0; i < M;)
            {
                NodeType v = curand(&localState) % u;
                double r = curand_uniform(&localState);

                if (r < p)
                {
                    //F_i(u) <- v
                    if (V.insertItem(i, v))
                    {
                        i++; //Success, next outgoing edge
                    }
                }
                else
                {
                    //F_i(u) <- F_k(v)
                    EdgeIndexType k = curand(&localState) % M;

                    NodeType F_k_v = G[v * M + k];

                    //F_k(v) is Known, Try to insert or retry
                    if (V.insertItem(i, F_k_v))
                    {
                        i++; //Success, next outgoing edge
                    }
                }
            }
        }
        free(V._leftChild);
        free(V._rightChild);
        state[id] = localState; //Restore Random Number State to the global memory
    }
}
