/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "cuPPA.cuh"
#include <helper_cuda.h>

cuPPA::cuPPA(NodeType nNodes, EdgeIndexType mEdgesPerNode, double p)
{
    _N = nNodes;
    _M = mEdgesPerNode;
    _p = p;

    _dGraph = NULL;
}

cuPPA::~cuPPA()
{
}

void cuPPA::initialize()
{
//    checkCudaErrors(cudaMalloc(&_dGraph, _N * _M * sizeof(NodeType)));
}

void cuPPA::generateGraphOld(ExperimentConfiguration conf)
{
    //void test_cuPPA(NodeType N, EdgeIndexType M, double p)
    //{
    //    cudaDeviceReset();
    //
    //    curandStateMRG32k3a *devMRGStates;
    //    checkCudaErrors(cudaMalloc((void ** ) &devMRGStates, 64 * 64 * sizeof(curandStateMRG32k3a)));
    //    unsigned long long seed = 1234;
    //    setup_curand_kernel<<<64, 64>>>(devMRGStates, seed);    //Creates 64*64 Streams!
    //
    //    int nThreads = THREADS_PER_BLOCK;
    //    int nBlocks = 2;
    //
    //    NodeType *dMem;
    //    EdgeIndexType *children;
    //    cuBSTVertex *vertices;
    //
    //    checkCudaErrors(cudaMalloc((void** ) &dMem, N * M * sizeof(NodeType)));
    //    checkCudaErrors(cudaMalloc((void** ) &children, 2 * N * M * sizeof(EdgeIndexType)));
    //    checkCudaErrors(cudaMalloc((void** ) &vertices, N * sizeof(cuBSTVertex)));
    //
    //    initialize_vertex_structures<<<nBlocks, nThreads>>>(N, M, dMem, children, vertices);
    //    initial_graph<<<nBlocks, nThreads>>>(M, vertices);
    ////    generate_graph_race<<<nBlocks, nThreads>>>(N, M, p, vertices, dMem, devMRGStates);
    //    generate_graph_avoid<<<nBlocks, nThreads>>>(N, M, p, vertices, dMem, devMRGStates);
    //
    //    cudaDeviceSynchronize();
    //
    //    NodeType *G = (NodeType*) malloc(N * M * sizeof(NodeType));
    //
    //    cudaMemcpy(G, dMem, N * M * sizeof(NodeType), cudaMemcpyDeviceToHost);
    //
    //    int* deg = (int*) calloc(N, sizeof(int));
    //
    //    for (int i = M; i < N; i++)
    //    {
    //        deg[i] = M;
    //        for (int j = 0; j < M; j++)
    //        {
    //            NodeType u = G[i * M + j];
    //
    //            if (u < 0)
    //            {
    //                printf("Problem: %d, %lld\n", i, u);
    //            }
    //            deg[u]++;
    //        }
    //    }
    //
    //    int* DD = (int*) calloc(N, sizeof(int));
    //
    //    for (int i = 0; i < N; i++)
    //    {
    //        if (deg[i] > 0)
    //        {
    //            DD[deg[i]]++;
    //        }
    //    }
    //
    //    printf("--------------------------------------\n");
    //    for (int i = 0; i < N; i++)
    //    {
    //        if (DD[i] > 0)
    //        {
    //            printf("%d\t%d\n", i, DD[i]);
    //        }
    //    }
    //    printf("--------------------------------------\n");
    //
    //    std::cout << "Ended" << std::endl;
    //    cudaDeviceSynchronize();
    //
    //}
}

void cuPPA::printDegreeDistribution(ExperimentConfiguration conf)
{
    if (conf._FLAG_SHOW_DEGREE_DISTRIBUTION)
    {
        printf("Generating Degree Distribution");

        NodeType *hG = (NodeType*) malloc(_N * _M * sizeof(NodeType));

        cudaMemcpy(hG, _dGraph, _N * _M * sizeof(NodeType), cudaMemcpyDeviceToHost);

        EdgeType* deg = (EdgeType*) calloc(_N, sizeof(EdgeType));

        for (NodeType uu = _M; uu < _N; uu++)
        {
            deg[uu] = _M; //Outgoing Edges
            for (int j = 0; j < _M; j++)
            {
                NodeType v = hG[uu * _M + j];

                if (v < 0)
                {
                    printf("Problem: %d, %d\n", uu, v);
                }
                else
                {
                    deg[v]++;
                }
            }
        }

        for (NodeType uu = 0; uu < _M; uu++)
        {
            cout << uu << "\t" << deg[uu] << endl;
        }

        EdgeType* DD = (EdgeType*) calloc(_N, sizeof(EdgeType));

        for (NodeType uu = 0; uu < _N; uu++)
        {
            if (deg[uu] > 0)
            {
                DD[deg[uu]]++;
            }
        }

        std::stringstream sstm;
        sstm << conf._outputDir << "/" << conf._distfile1;
        std::ofstream ddFile(sstm.str().c_str());

        ddFile << "Degree" << "," << "Count" << endl;
        for (NodeType uu = 0; uu < _N; uu++)
        {
            if (DD[uu] > 0)
            {
                ddFile << uu << "," << DD[uu] << endl;
            }
        }
        ddFile.close();

        free(deg);
        free(DD);
    }
}
