/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "Tests1.cuh"
#include "../Graph/Vertex/CPU/BSTVertex.h"
#include "../Graph/Vertex/CUDA/cuBSTVertex.cuh"
#include "../Common/Random.hpp"
#include "../Common/cuRandom.cuh"
#include <helper_cuda.h>
#include "../PA/SPA/SerialCopyModel.h"
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void kernel_atomic_add()
{

    int id = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ int sharedVariable;
    sharedVariable = 0;
    __syncthreads();

//    printf("Thread Id: %d\n", id);

    printf("%d\n", sharedVariable);

    if (id / 8)
    {
        int k = atomicAdd(&sharedVariable, 1);
    }

    __syncthreads();

    printf("Thread Id: %d\t%d\n", id, sharedVariable);

}

void testAtomicAdd()
{
    cudaDeviceReset();
    kernel_atomic_add<<<2,8>>>();
    cudaDeviceSynchronize();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void testVertexStructureInCUDAKernel(curandStateMRG32k3a *state, NodeType N, EdgeIndexType M, double p, NodeType* A)
{
    cuBSTVertex V(M);

    int loop_count = 0;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;

    curandStateMRG32k3a localState = state[tid]; //Save Random Number State to the global memory

    double B = (double) N / (double) blockDim.x;

    NodeType lb = tid * B;
    lb = lb < M + 1 ? M + 1 : lb;

    NodeType ub = (tid + 1) * B - 1;
    ub = ub > N - 1 ? N - 1 : ub;

    printf("%d: %f, %ld, %ld\n", tid, B, lb, ub);

    NodeType u = 17;

//    NodeType uu = u;

    int qIndex = 0;

    for (NodeType uu = lb; uu < ub; uu++)
//    for (NodeType k = 0; k < 10; k++)
    {
        printf("O: %d: %ld\n", tid, uu);

        //Step 1: Computating the
        for (EdgeIndexType i = 0; i < M;)
        {
            NodeType v = curand(&localState) % uu;
            loop_count++;

            double r = curand_uniform(&localState);

            if (r < p)
            {
                if (V.insertItem(i, v) >= 0)
                    i++;
            }
            else
            {
//                int currentLocation = atomicAdd(&qIndex, 1);

                //Try to resolve the items first.
                //If not possible, put it into the waiting list.
            }
        }
//        printf("L:%d\t%ld\t%d\n", tid, u, loop_count);

        loop_count = 0;
        V.clear();
    }
    state[tid] = localState; //Restore Random Number State to the global memory

    free(V._A);
    free(V._rightChild);
    free(V._leftChild);

}

void testVertexStructureInCUDA(EdgeIndexType M)
{
    cudaDeviceReset();

    curandStateMRG32k3a *devMRGStates;
    checkCudaErrors(cudaMalloc((void ** ) &devMRGStates, 64 * 64 * sizeof(curandStateMRG32k3a)));
    unsigned long long seed = 1234;
    setup_curand_kernel_MRG32k3a<<<64, 64>>>(devMRGStates, seed);    //Creates 64*64 Streams!

    int nThreads = 16;

    NodeType N = 1024;

    NodeType *dMem;
//    checkCudaErrors(cudaMalloc((void** ) &dMem, N * M * sizeof(NodeType)));

//    getLastCudaError("testVertexStructureInCUDAKernel() execution FAILED\n");

//    NodeType *hMem = (NodeType*) calloc(N * M, sizeof(NodeType));

//    checkCudaErrors(cudaMemcpy(dMem, hMem, N * M * sizeof(NodeType), cudaMemcpyHostToDevice));

//    testVertexStructureInCUDAKernel<<<1,nThreads>>>(devMRGStates, N, M, 0.5, dMem );

//    checkCudaErrors(cudaMemcpy(hMem, dMem, N * M * sizeof(NodeType), cudaMemcpyDeviceToHost));

//    for (int i = 0; i < nThreads; i++)
//    {
//        std::cout << hMem[i] << std::endl;
//    }

    std::cout << "Ended" << std::endl;
    cudaDeviceSynchronize();

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void testVertexStructure(EdgeIndexType M)
{
    NodeType U = M + 1 + 16;
    long long total = 0;

    Timer t;
    t.start();

    for (NodeType u = M + 1; u < U; u++)
    {
        BSTVertex V(M);
//        HashVertex V(M);

        int loop_count = 0;

        for (EdgeIndexType i = 0; i < M;)
        {
            NodeType v = Random::Uniform() % u;

            if (V.insertItem(v) >= 0)
                i++;

            loop_count++;
        }
//        std::cout << "Loop Count:\t";
        std::cout << loop_count << endl;
        total += loop_count;
    }
    t.stop();

    std::cout << total << "\t" << t.getElapsedTimeInSec() << endl;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void test_SerialCopyModel(NodeType N, EdgeIndexType M, double p)
{
    SerialCopyModel scm(N, M, p);
    ExperimentConfiguration conf;
    conf._n1 = N;
    conf._x1 = M;
    conf._alpha = p;
    conf._outputDir = "./";
    conf._distfile1 = "dd.csv";
    scm.generateGraph(conf);
    printf("Ended\n");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void kernel_FindMaxGlobally(int data, int X[])
{
    __shared__ EdgeType WQ[6 * 1024];

//    WQ[threadIdx.x] = blockDim.x * gridDim.x - threadIdx.x - blockIdx.x * blockDim.x;
    WQ[threadIdx.x] = threadIdx.x + blockIdx.x * blockDim.x;
    __syncthreads();

    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (threadIdx.x < s)
        {
            WQ[threadIdx.x] = max(WQ[threadIdx.x], WQ[threadIdx.x + s]);
        }
        __syncthreads();
    }

    X[blockIdx.x] = WQ[0];
}

__global__ void kernel_FindMAX(int X[])
{
    __shared__ EdgeType WQ[6 * 1024];

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

    X[0] = WQ[0];
}

void testFindMaxGlobally()
{
    int N_BLOCK1 = 32;
    int N_BLOCK = 13;
    int* dPtr = NULL;

    printf("TEST\n");

    cout << dPtr << endl;
    int data = 10;

    cudaDeviceReset();

    int *hPtr = (int*) malloc(N_BLOCK1 * sizeof(int));

    for (int i = 0; i < N_BLOCK1; i++)
    {
        hPtr[i] = 100 + i;
    }

    checkCudaErrors(cudaMalloc(&dPtr, N_BLOCK1 * sizeof(int)));
    checkCudaErrors(cudaMemset(dPtr, 1, N_BLOCK1 * sizeof(int)));

    checkCudaErrors(cudaMemcpy(dPtr, hPtr, N_BLOCK1 * sizeof(int), cudaMemcpyHostToDevice));

    kernel_FindMaxGlobally<<<N_BLOCK, 512>>>(data, dPtr);
    kernel_FindMAX<<<1, N_BLOCK1>>>(dPtr);

    checkCudaErrors(cudaMemcpy(hPtr, dPtr, N_BLOCK1 * sizeof(int), cudaMemcpyDeviceToHost));

    cout << "\t" << endl;

    for (int i = 0; i < 1; i++)
    {
        cout << "Block:  " << i << "\t" << hPtr[i] << endl;
    }

    cudaDeviceSynchronize();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
