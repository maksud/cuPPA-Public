/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "SerialBA.h"
#include <iostream>
#include <iterator>
#include <algorithm>
#include "../../Common/Random.hpp"
#include <string>
#include <sstream>
#include "PAHelper.hpp"
#include "../../Common/Timer.h"
#include <set>

SerialBA::SerialBA(LongLong n, LongLong m)
{
    _N = n;
    _M = m;
}

NodeType* SerialBA::initializeNetwork(std::ofstream& log)
{
#if SerialBA_LOG
    log << "Required Memory for Graph: " << (sizeof(NodeType) * 2 * _N * _M) / (1024.0 * 1024.0) << std::endl;
#endif
    size_t size = 2 * (uint64_t) _N * _M * sizeof(NodeType);
    NodeType* G = (NodeType*) malloc(size);
    if (G == NULL)
    {
        log << "Not Enough Memory!" << endl;
        return NULL;
    }
    log << "Memory Allocated" << endl;
    return G;
}

SerialBA::~SerialBA()
{
}

void SerialBA::clearNetwork(NodeType *G)
{
    free(G);
}

int SerialBA::findNode(const NodeType degree[], const NodeType N, NodeType R)
{
    for (NodeType i = 0; i < N; i++)
    {
        R -= degree[i];
        if (R <= 0)
            return (i);
    }
    return (-1);
}

inline void SerialBA::randomSet(NodeType nodes[], std::set<NodeType>& target, EdgeType Z)
{
    target.clear();
    while (target.size() < _M)
    {
        target.insert(nodes[Random::UniformMax(Z)]);
    }
}

void SerialBA::generateGraphNetworkX(ExperimentConfiguration conf)
{
    string METHOD = "SBA";

    std::cout << METHOD << std::endl;

    //LOG
    std::stringstream sstm;
    sstm << conf._outputDir << METHOD << ".log";
    std::ofstream log(sstm.str().c_str());

    log << METHOD << std::endl;
    log << _N << std::endl;
    log << _M << std::endl;

    /*****************************************************************************/
    NodeType* _G = initializeNetwork(log);

    Timer tComputation;

    tComputation.start();

    for (NodeType i = 0; i < _M; i++)
    {
        for (EdgeIndexType j = 0; j < _M; j++)
        {
            _G[i * _M * 2 + j] = i;
            _G[i * _M * 2 + _M + j] = i;
        }
    }

#if 0
    set<NodeType> target;
    for (NodeType i = 0; i < _M; i++)
    {
        repeated_nodes[i * 2] = i;
        repeated_nodes[i * 2 + 1] = i;
        target.insert(i);
    }
#else
    BSTVertex uniqueNodes(_M);
    for (NodeType i = 0; i < _M; i++)
    {
        uniqueNodes.insertItem(i, i);
    }
#endif

    //Working Loop
    EdgeType Z = 2 * _M * _M;
    for (NodeType i = _M; i < _N; i++)
    {
#if SerialBA_DEBUG
        std::cout << "\n";
        std::cout << "Processing Node: " << (i+1) << std::endl;
        PAHelper::print_array("repeated_nodes: ", repeated_nodes, Z);
        PAHelper::print_array("target: ", target, _M);
#endif

#if 0
        set<NodeType>::iterator it = target.begin();
        for (NodeType j = 0; j < _M; j++)
        {
            repeated_nodes[Z + j] = *it;
            repeated_nodes[Z + _M + j] = i;
            it++;
        }
        Z += 2 * _M;

        randomSet(repeated_nodes, target, Z);
#else
        for (NodeType j = 0; j < _M; j++)
        {
//            _G[Z + j] = uniqueNodes.Adj(j);
//            _G[Z + _M + j] = i;

            _G[Z + 2 * j] = uniqueNodes.Adj(j);
            _G[Z + 2 * j + 1] = i;
        }
        Z += 2 * _M;

        //BST Based Approach?
        uniqueNodes.clear();

        for (NodeType j = 0; j < _M;)
        {
            if (uniqueNodes.insertItem(j, _G[Random::UniformMax(Z)]))
            {
                j++;
            }
        }

#endif
    }
    tComputation.stop();
    double computationTime = tComputation.getElapsedTimeInSec();

#if SerialBA_LOG
    std::cout << "Computation Time : " << computationTime << std::endl;
    log << "Computation Time : " << computationTime << std::endl;
//	Utilities::logMemUsage(&log, "");
#else
    std::cout << "VM Used          : " << PAHelper::getUsedVirtualMemory() << std::endl;
    std::cout << "Phys Used        : " << PAHelper::getUsedPhysicalMemory() << std::endl;
#endif

#if SerialBA_OUTPUT_FILE
    gettimeofday(&start, 0);
    PAHelper::output_ba_graph1(_G, _N, _M, 0, 0, outputDir);
    gettimeofday(&end, 0);
    long outputTime = Utilities::difftime(end, start);
    log << "Output Time      : " << outputTime << std::endl;
#endif

    if (conf._FLAG_SHOW_DEGREE_DISTRIBUTION)
    {
        PAHelper::outputSBADegreeDistribution(_G, _N, _M, conf._outputDir, conf._distfile1.c_str(), true);
    }

    clearNetwork(_G);
}

void SerialBA::generateGraph(ExperimentConfiguration conf)
{
    generateGraphNetworkX(conf);
}

