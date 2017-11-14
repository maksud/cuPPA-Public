/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "main.cuh"
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>

#include "Graph/Vertex/CPU/BSTVertex.h"
#include "Graph/Vertex/CUDA/cuBSTVertex.cuh"
#include "Graph/Vertex/CPU/HashVertex.h"
#include "Graph/Vertex/CPU/LinearVertex.h"
#include "Common/Random.hpp"
#include "Common/cuRandom.cuh"
#include "Common/Timer.h"
#include <cstdio>
#include <helper_cuda.h>

#include "PA/cuPPA/cuPPA.cuh"
#include "PA/SPA/SerialCopyModel.h"
#include "PA/SPA/SerialBA.h"

#include "Common/CommonHelper.cuh"

#include "Utility/CLI/CLI.hpp"
#include "Utility/ExperimentConfiguration.hpp"

#include "Test/Tests1.cuh"

/**
 * Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 */
static void CheckCudaErrorAux(const char *file, unsigned line, const char *statement, cudaError_t err)
{
    if (err == cudaSuccess)
        return;
    std::cerr << statement << " returned " << cudaGetErrorString(err) << "(" << err << ") at " << file << ":" << line << std::endl;
    exit(1);
}

void syntheticNetworkGeneration(ExperimentConfiguration conf, int argc, char *argv[])
{
    int rank = 0, P = 1;

    if (rank == 0)
    {
        printf("Method:\t%s\t", conf._algorithm.c_str());
        printf("Load Balancing:\t%s\t", conf._loadBalancing.c_str());
        printf("Weight:\t%s\n", conf._weight1.c_str());
        printf("N:\t%ld\t", conf._n1);
        printf("M:\t%u\n", conf._x1);
        printf("Min Degree:\t%d\t", conf._d1min);
        printf("Max Degree:\t%d\t", conf._d1max);
        printf("Avg Degree:\t%d\n", conf._d1avg);
        printf("Output:\t%s\n", conf._outputDir.c_str());
        printf("Probability:\t%f\t", conf._probability);
        printf("Processors:\t%d\t", P);
        printf("Gamma:\t%f\t", conf._gamma);
        printf("a:\t%f\n", conf._alpha);
        printf("MSG Buffer:\t%ld\t", conf._msgBufferSize);
        printf("MPI Buffer:\t%ld\t", conf._mpiBufferSize);
        printf("LB Factor:\t%f\t", conf._lbFactor);
        printf("Num Stages:\t%ld\t", conf._stages);
        printf("Write To File:\t%d\n", conf._write);

//      printf("Virtual  Memory Used: %d KB\n", Utilities::getUsedVirtualMemory());
//      printf("Physical Memory Used: %d KB\n", Utilities::getUsedPhysicalMemory());

    }

    GenerationAlgorithm method;

//    LoadType weight1;
//    LoadType weight2;
//    LoadBalancingScheme lbScheme;
//    WeightType* iW1 = NULL;
//    WeightType* iW2 = NULL;

    Methods methods;
    method = methods.ParseAlgorithm(conf._algorithm);

//    weight1 = methods.ParseLoadType(conf._weight1);
//    weight2 = methods.ParseLoadType(conf._weight2);
//    lbScheme = methods.ParseLoadBalancingScheme(conf._loadBalancing);

    switch (method)
    {
    case METHOD_SERIAL_COPY_MODEL:
    {
        SerialCopyModel ba(conf._n1, conf._x1, conf._alpha);
        ba.generateGraph(conf);
        break;
    }
    case METHOD_SERIAL_BA_NETWORKX:
    {
        SerialBA ba(conf._n1, conf._x1);
        ba.generateGraph(conf);
        break;
    }
    case METHOD_CUPPA:
    {
        cuPPA cuppa(conf._n1, conf._x1, conf._alpha);
        cuppa.generateGraph(conf);
        break;
    }
    case METHOD_CUPPA_ADAPTIVE:
    {
        cuPPA cuppa(conf._n1, conf._x1, conf._alpha);
        cuppa.generateGraphAdaptive(conf);
        break;
    }
    default:
        throw std::logic_error(__FILE__ ": No method specified.");
    }

}

ExperimentConfiguration parseExperiments(int argc, char *argv[])
{
    ExperimentConfiguration expt;

    CLI::App app("cuPSyNet");

    {
        expt._n1 = 1000;
        expt._x1 = 4;
        expt._d1min = 0;
        expt._d1max = 0;
        expt._d1avg = 0;

        expt._weight1 = "none";
        expt._distfile1 = "degreeDistribution.txt";

        expt._avgcc1 = "";
        expt._probabilityMatrixFile1 = "";

        expt._alpha = 0.5;
        expt._gamma = 2.7;
    }
    {
        expt._n2 = 1000;
        expt._x2 = 4;
        expt._d2min = 0;
        expt._d2max = 0;
        expt._d2avg = 0;

        expt._weight2 = "none";
        expt._distfile2 = "degreeDistribution.txt";

        expt._avgcc2 = "";
        expt._probabilityMatrixFile2 = "";
    }
    {
        expt._outputDir = "./";
        expt._msgBufferSize = 0;
        expt._mpiBufferSize = 0;
    }
    {
        expt._lbFactor = 1;
        expt._probability = 0.1;
        expt._stages = 1000;
        expt._loadBalancing = "none";
        expt._write = false;
    }
    {
        expt._blockSize = 2;
        expt._threadsPerBlock = 32;
        expt._warpSize = 32;
    }
    {
        expt._msgBufferSize = 1024;
        expt._mpiBufferSize = 1024;
    }

    std::string conf_ini;

    app.add_option("-c,--config-file", conf_ini, "Configuration Filename");
    {
        app.add_option("-A,--algorithm", expt._algorithm, "Algorithm");
        app.add_option("-l,--loadbalancing", expt._loadBalancing, "Load Balancing");
    }
    {
        app.add_option("-n,--n1", expt._n1, "number of vertices");
        app.add_option("-x,--x1", expt._x1, "x1");
        app.add_option("--weight1", expt._weight1, "Weight 1");
    }
    {
        app.add_option("-N,--n2", expt._n2, "number of vertices 2");
        app.add_option("-X,--x2", expt._x2, "x2");
        app.add_option("--weight2", expt._weight2, "Weight 2");
    }
    {
        app.add_option("-a,--alpha", expt._alpha, "Probability of creating a direct edge");
    }
    {
        app.add_option("-s,--stages", expt._stages, "#of Stages");
    }
    {
        //MPI-Related
        app.add_option("--msg-buffer-size", expt._msgBufferSize, "Buffer Size");
        app.add_option("--mpi-buffer-size", expt._mpiBufferSize, "MPI Buffer Size");
    }
    {
        //CUDA Related
        app.add_option("-b,--blocks", expt._blockSize, "# of Blocks");
        app.add_option("-t,--threads", expt._threadsPerBlock, "# of Threads per Block");
        app.add_option("-w,--warps", expt._warpSize, "Warp Size");
    }
    {
        //Degree Distribution
        app.add_option("--show-dd", expt._FLAG_SHOW_DEGREE_DISTRIBUTION, "Export degree distribution?");
        app.add_option("--dd-filename", expt._distfile1, "Degree Distribution Filename");
        app.add_option("--show-time-n-stages", expt._N_STAGE_TIMINGS, "Number of step timings to show");
    }
    try
    {
        app.parse(argc, argv);
    } catch (const CLI::Error &e)
    {
        cout << "Exception" << e.what() << endl;
        exit(0);
    }

    if (app.count("--config-file"))
    {
//        cout << "Loading Config File: " << conf_ini << endl;
        expt.load(conf_ini.c_str());
    }

    return expt;
}

void printWelcomecuPPA()
{
    cout << "*********************************************************************" << endl;
    cout << endl;
    cout << "               ______  ______    ___  " << endl;
    cout << "               | ___ \\ | ___ \\  / _ \\ " << endl;
    cout << "  ___   _   _  | |_/ / | |_/ / / /_\\ \\" << endl;
    cout << " / __| | | | | |  __/  |  __/  |  _  |" << endl;
    cout << "| (__  | |_| | | |     | |     | | | |" << endl;
    cout << " \\___|  \\__,_| \\_|     \\_|     \\_| |_/" << endl;
    cout << endl;
    cout << "*********************************************************************" << endl;
    cout << endl;
}

int main(int argc, char *argv[])
{
    printWelcomecuPPA();

    ExperimentConfiguration expt = parseExperiments(argc, argv);
    syntheticNetworkGeneration(expt, argc, argv);
    return 0;
}

