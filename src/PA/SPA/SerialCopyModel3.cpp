/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "SerialCopyModel.h"
#include "../../Common/Random.hpp"
#include "PAHelper.hpp"
#include "../../Graph/Vertex/CPU/BSTVertex.h"

void SerialCopyModel::generateGraphCopyModelBST(ExperimentConfiguration conf)
{
    Timer t;
    string METHOD = "SCM";

    std::cout << METHOD << std::endl;

    //LOG
    std::stringstream sstm;
    sstm << conf._outputDir << METHOD << ".log";
    std::ofstream log(sstm.str().c_str());

    log << METHOD << std::endl;
    log << _N << std::endl;
    log << _M << std::endl;

    t.start();
    /*****************************************************************************/
    NodeType *_G = initializeNetwork(log, -1);
    BSTVertex uniqueNodes(_M);

//    BSTVertex bb(_M);

#if SerialCopyModel_RETRY_COUNT
    int retry_count = 0;
#endif

    for (EdgeIndexType i = 0; i < _M; i++)
    {
        for (EdgeIndexType j = 0; j < _M; j++)
        {
            _G[i * _M + j] = i;
        }
    }

    for (EdgeIndexType j = 0; j < _M; j++)
    {
        _G[_M * _M + j] = j;
    }

    // From M+1 nodes apply the Copy Model
    for (EdgeType i = _M + 1; i < _N; i++)
    {
        uint64_t Z = _M * i; //Total Edges upto this node

        uniqueNodes.clear();
        for (EdgeIndexType j = 0; j < _M; j++)
        {
            NodeType target = -1;
            do
            {
                EdgeType R = Random::UniformMax(Z);
                NodeType node_index = R / _M;
                EdgeIndexType edge_index = R % _M;
                //
                double coin = Random::Uniform01();
                if (coin < _p)
                {
                    target = node_index;
                }
                else
                {
                    target = _G[R];
                }
            } while (!uniqueNodes.insertItem(j, target));
            _G[Z + j] = uniqueNodes[j];
        }
    }
    t.stop();
    double computationTime = t.getElapsedTimeInSec();

    t.reset();
    t.start();
#if SerialCopyModel_OUTPUT_FILE
    std::stringstream sstmGraphFile;
    sstmGraphFile << conf._outputDir << METHOD << ".edges";
    PAHelper::writeGraphToDisk(_G, _N, _M, sstmGraphFile.str());
#endif
    t.stop();
    double outputTime = t.getElapsedTimeInSec();

#if SerialCopyModel_LOG
    std::cout << "Computation Time : " << computationTime << std::endl;
    log << "Computation Time : " << computationTime << std::endl;
    log << "Output Time      : " << outputTime << std::endl;
    Utilities::logMemUsage(&log, "");
#else
    std::cout << "Output Time      : " << outputTime << std::endl;
    std::cout << "VM Used          : " << PAHelper::getUsedVirtualMemory() << std::endl;
    std::cout << "Phys Used        : " << PAHelper::getUsedPhysicalMemory() << std::endl;
#endif

    if (conf._FLAG_SHOW_DEGREE_DISTRIBUTION)
    {
        PAHelper::outputSCMDegreeDistribution(_G, _N, _M, conf._outputDir, conf._distfile1.c_str());
    }
}
