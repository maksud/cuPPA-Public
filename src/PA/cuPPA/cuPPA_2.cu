/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "cuPPA.cuh"
#include <helper_cuda.h>

void cuPPA::generateGraphAdaptive(ExperimentConfiguration conf)
{
    cudaDeviceReset();

    _N = conf._n1;
    _M = conf._x1;
    _p = conf._alpha;

#if WAITING_QUEUE_IN_GLOBAL_MEMORY
    const int WQ_MAX_SIZE = 1024;
#else
    const int WQ_MAX_SIZE = 48;
#endif
//    const int nBlocks = getSMCount();
    const int nBlocks = conf._blockSize;
    const int MAX_THREADS_PER_BLOCK = conf._threadsPerBlock;
    const int MAX_QUEUE_CAPACITY = WQ_MAX_SIZE * 1024 / sizeof(EdgeType);
    const NodeType MAX_VERTICES_IN_BST = 32;

    //INITIAL SETTINGS
    int nThreads = min(Utilities::lowerPower2(MAX_QUEUE_CAPACITY / _M), MAX_THREADS_PER_BLOCK);
    float CURRENT_QUEUE_SIZE = MAX_QUEUE_CAPACITY / nThreads;
    int VERTICES_PER_THREAD = MAX_QUEUE_CAPACITY / nThreads / _M;
    int REQUIRED_QUEUE = _M * VERTICES_PER_THREAD;
    const int MAX_BLOCK_SIZE = min(MAX_THREADS_PER_BLOCK * nBlocks * MAX_VERTICES_IN_BST, _N - _M - 1);

#if USE_DEFAULT_CURAND
    curandState *devStates;
    checkCudaErrors(cudaMalloc((void ** ) &devStates, nBlocks * MAX_THREADS_PER_BLOCK * sizeof(curandState)));
    //    unsigned long long seed = 1234;
    unsigned long long seed = std::time(0);
    setup_curand_kernel<<<nBlocks, MAX_THREADS_PER_BLOCK>>>(devStates, seed); //Creates nBlocks * MAX_THREADS_PER_BLOCK Streams!
#else
    curandStateMRG32k3a *devStates;
    checkCudaErrors(cudaMalloc((void ** ) &devStates, nBlocks * MAX_THREADS_PER_BLOCK * sizeof(curandStateMRG32k3a)));
    unsigned long long seed = 1234;
    setup_curand_kernel_MRG32k3a<<<nBlocks, MAX_THREADS_PER_BLOCK>>>(devStates, seed);    //Creates nBlocks * nThreads Streams!
#endif

    checkCudaErrors(cudaMalloc((void** ) &_dGraph, _N * _M * sizeof(NodeType))); //Whole Graph Size

    cuBSTVertex *F;
    checkCudaErrors(cudaMalloc((void** ) &F, MAX_BLOCK_SIZE * sizeof(cuBSTVertex))); //cuBSTVertex

    EdgeIndexType * children;
    checkCudaErrors(cudaMalloc((void** ) &children, 2 * MAX_BLOCK_SIZE * _M * sizeof(EdgeIndexType))); //Left and Right Children

    int *dMaxQueueSizes;
    int hMaxQueueSize = 0;

#if WAITING_QUEUE_IN_GLOBAL_MEMORY
    EdgeType *dQ; //Shared Memory Size
    checkCudaErrors(cudaMalloc((void** ) &dQ, WQ_MAX_SIZE * 1024 * nBlocks));//Waiting Queue
    checkCudaErrors(cudaMalloc((void** ) &dMaxQueueSizes, nThreads * nBlocks * sizeof(int)));
    cudaMemset(dMaxQueueSizes, 0, nThreads * nBlocks * sizeof(int));
#else
    checkCudaErrors(cudaMalloc((void** ) &dMaxQueueSizes, 32 * sizeof(int)));
    cudaMemset(dMaxQueueSizes, 0, 32 * sizeof(int));
#endif

#if SHOW_QUEUE_OUTPUT
    std::ofstream maxQFile("Q-Adaptive.csv");
    maxQFile << "Round,MaxQ\n";
#endif

#if SHOW_INTERNAL_TIMINGS
    std::ofstream stepwiseTimingFile("StepwiseTiming-Adaptive.csv");
    stepwiseTimingFile << "Round,nBlocks,nThreads,Edges,MaxQ,Initialization Time,Computation Time\n";
#endif

    printf("--------------------------------------------------------------------------------------------------------------------------------------------\n");
    cout << "BLOCKS             : \t" << nBlocks << "\t| " << "MAX_THREADS_PER_BLOCK: \t" << MAX_THREADS_PER_BLOCK << "\t| " << "nThreads     : " << nThreads << endl;
    printf("--------------------------------------------------------------------------------------------------------------------------------------------\n");
    cout << "MAX_QUEUE_CAPACITY : \t" << MAX_QUEUE_CAPACITY << "\t| " << "REQUIRED_QUEUE       : \t" << REQUIRED_QUEUE << "\t| " << "CURRENT_QUEUE: " << MAX_QUEUE_CAPACITY / nThreads << endl;
    printf("--------------------------------------------------------------------------------------------------------------------------------------------\n");
    cout << "VERTICES_PER_THREAD: \t" << VERTICES_PER_THREAD << "\t| " << "MAX_BLOCK_SIZE       : \t" << MAX_BLOCK_SIZE << "\t| " << "STEPS        : \t" << "N/A" << endl;
    printf("--------------------------------------------------------------------------------------------------------------------------------------------\n");
    cout << "Random: \t" << (nBlocks * MAX_THREADS_PER_BLOCK * sizeof(curandState) / 1024.0 / 1024.0) << "MB" << "\t| " //
            << "Graph: \t" << (_N * _M * sizeof(NodeType) / 1024.0 / 1024.0) << "MB" << "\t| " //
            << "BST: \t" << (MAX_BLOCK_SIZE * sizeof(cuBSTVertex) / 1024.0 / 1024.0) << "MB" << "\t| " //
            << "Children: \t" << (2 * MAX_BLOCK_SIZE * _M * sizeof(EdgeIndexType) / 1024.0 / 1024.0) << "MB" << "\t| " //
            << endl;
    printf("--------------------------------------------------------------------------------------------------------------------------------------------\n");

    Timer t1;
    t1.start();
    initial_graph_rrp<<<nBlocks, nThreads>>>(_M, _dGraph);

    int iteration_index = -1;

    NodeType u = _M + 1;
    int BLOCK_SIZE = VERTICES_PER_THREAD * nThreads * nBlocks;
    bool FINAL_TUNING = false;

    while (u < _N)
    {

#if SHOW_INTERNAL_TIMINGS & 0
        cout << endl << "Blocks:\t" << nBlocks << "\tThreads:\t" << nThreads << endl << endl;
        cout << "Block Size: " << BLOCK_SIZE << endl;
        cout << "Current Queue Size: " << CURRENT_QUEUE_SIZE << endl;
#endif

        iteration_index++;
        NodeType v = u + BLOCK_SIZE - 1;
        v = v > _N - 1 ? _N - 1 : v;
        NodeType B = v - u + 1;

#if SHOW_INTERNAL_TIMINGS & 0
        if (iteration_index < conf._N_STAGE_TIMINGS)
        {
            printf("%d\tProcessing: %lld\t%lld\tTotal: %lld\n", iteration_index, u, v, B);
        }
#endif
        if (u < _N)
        {
#if SHOW_INTERNAL_TIMINGS
            Timer tVS;
            tVS.start();
#endif

#if 0
            //NOT SUPER FAST as expected
            //Vertex Setup can be done super fast without any dependency!
            int nThreadsVertexSetup = 512;
            int nBlocksVertexSetup = (B + nThreadsVertexSetup - 1) / nThreadsVertexSetup;
            initialize_vertex_structures_rrp_2<<<nBlocksVertexSetup, nThreadsVertexSetup>>>(u, B, M, _dGraph, children, F);
#else
#if 0
            //NOT FASTER, but expected it to be.
            cudaMemset(children, 0, 2 * B * _M * sizeof(EdgeIndexType));
            cudaMemset(_dGraph + (u * M), 0, B * _M * sizeof(NodeType));
            initialize_vertex_structures_rrp_fast<<<nBlocks, nThreads>>>(u, B, _M, _dGraph, children, F);
#else
            //Currently the fastest implementation
            initialize_vertex_structures_rrp<<<nBlocks, nThreads>>>(u, B, _M, _dGraph, children, F);
#endif
#endif

#if SHOW_INTERNAL_TIMINGS
            cudaDeviceSynchronize();
            tVS.stop();
#endif

#if SHOW_INTERNAL_TIMINGS
            Timer tGn;
            tGn.start();
#endif

#if WAITING_QUEUE_IN_GLOBAL_MEMORY
            generate_graph_rrp_by_batch_execute_copy_model<<<nBlocks, nThreads>>>(u, B, _M, _p, F, _dGraph, devStates, dQ, dMaxQueueSizes);
            generate_graph_rrp_by_batch_resolve_queue<<<nBlocks, nThreads>>>(u, B, _M, _p, F, _dGraph, devStates, dQ, dMaxQueueSizes);
#else
            generate_graph_rrp_by_batch<<<nBlocks, nThreads>>>(u, B, _M, _p, F, _dGraph, devStates, dMaxQueueSizes);
            find_max_queue_size<<<1, 32>>>(dMaxQueueSizes);
#endif

#if SHOW_INTERNAL_TIMINGS
            cudaDeviceSynchronize();
            tGn.stop();
#endif

            cudaMemcpy(&hMaxQueueSize, dMaxQueueSizes, sizeof(int), cudaMemcpyDeviceToHost);
#if SHOW_QUEUE_OUTPUT
            if (iteration_index < conf._N_STAGE_TIMINGS)
            {
                int Cap = WQ_MAX_SIZE * 1024 / nThreads / sizeof(EdgeType);
                //printf("************************************************************************************                      MAX:\t%d/%d\n", hMaxQueueSize, Cap);
                maxQFile << (iteration_index + 1) << "," << hMaxQueueSize << endl;
            }
#endif

            if (!FINAL_TUNING && ((hMaxQueueSize == 0) | (CURRENT_QUEUE_SIZE / hMaxQueueSize > 2)))
            {
                if (nThreads < MAX_THREADS_PER_BLOCK)
                {
                    nThreads *= 2;
                    CURRENT_QUEUE_SIZE = CURRENT_QUEUE_SIZE / 2;
                }

#if SHOW_INTERNAL_TIMINGS & 0
                printf("********        TUNED PARAMETERS      *************** New BLOCK SIZE: %d\tNew Threads: %d\n", BLOCK_SIZE + NEXT_BLOCK_ADDITION, nThreads);
#endif

                u += BLOCK_SIZE;
                BLOCK_SIZE = min(BLOCK_SIZE * 2, MAX_BLOCK_SIZE); //Do NOT CROSS?

                if (nThreads == MAX_THREADS_PER_BLOCK && BLOCK_SIZE == MAX_BLOCK_SIZE)
                    FINAL_TUNING = true;
            }
            else
            {
                u += BLOCK_SIZE;
            }

#if SHOW_INTERNAL_TIMINGS
            if (iteration_index < conf._N_STAGE_TIMINGS)
            {
                //printf("Step:\t%d\tVertex Setup Time:\t%f (%3.2f Mps)\t\tGeneration Time:\t%f (%3.2f Mps)\n", iteration_index, tVS.getElapsedTimeInMilliSec(), B / tVS.getElapsedTimeInMilliSec() / 1024.0, tGn.getElapsedTimeInMilliSec(),
                //        B / tGn.getElapsedTimeInMilliSec() / 1024.0);
                stepwiseTimingFile << (iteration_index + 1) << "," << nBlocks << "," << nThreads << "," << B * _M << "," << hMaxQueueSize << "," << tVS.getElapsedTimeInSec() << "," << tGn.getElapsedTimeInSec() << endl;
            }
#endif
        }
    }

#if SHOW_QUEUE_OUTPUT
    maxQFile.close();
#endif

#if SHOW_INTERNAL_TIMINGS
    stepwiseTimingFile.close();
#endif

    cudaDeviceSynchronize();
    t1.stop();

    printf("Computation Time: %f\n", t1.getElapsedTimeInSec());
    printf("\n\n");

    if (conf._FLAG_SHOW_DEGREE_DISTRIBUTION)
    {
        printDegreeDistribution(conf);
    }

    cudaFree(devStates);
    cudaFree(_dGraph);
    cudaFree(children);
    cudaFree(F);

    std::cout << "Ended" << std::endl;
    cudaDeviceSynchronize();
}
