/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "SerialCopyModel.h"
#include <iostream>
#include <iterator>
#include <algorithm>
#include "../../Common/Random.hpp"
#include <string>
#include <sstream>
#include "../../Common/Timer.h"
#include <set>

/*
 * serialba.cpp
 *
 *  Created on: Feb 1, 2012
 *      Author: maksud
 */

SerialCopyModel::SerialCopyModel(NodeType n, EdgeIndexType m, double p)
{
    _N = n;
    _M = m;
    _p = p;
}

NodeType* SerialCopyModel::initializeNetwork(std::ofstream& log, int init_val = -1)
{
#if SerialCopyModel_LOG
    log << "Required Memory for Graph: " << (sizeof(NodeType) * _N * _M) / (1024.0 * 1024.0) << std::endl;
#endif
    NodeType* G = (NodeType*) malloc(_N * _M * sizeof(NodeType));
    if (G == NULL)
    {
        log << "Not Enough Memory!" << endl;
        return NULL;
    }
    std::fill_n(G, _M * _N, init_val);
    log << "Memory Allocated" << endl;
    return G;
}

SerialCopyModel::~SerialCopyModel()
{
}

void SerialCopyModel::clearNetwork(NodeType *G)
{
    free(G);
}

int SerialCopyModel::findNode(const NodeType degree[], const NodeType N, NodeType R)
{
    for (NodeType i = 0; i < N; i++)
    {
        R -= degree[i];
        if (R <= 0)
            return (i);
    }
    return (-1);
}

void SerialCopyModel::generateGraph(ExperimentConfiguration conf)
{
    generateGraphCopyModelBST(conf); //Faster for large M also memory efficient
}

