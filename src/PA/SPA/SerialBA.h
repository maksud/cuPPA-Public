/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef SERIALBA_ALGO_H_
#define SERIALBA_ALGO_H_

#include <string>
#include <set>
#include <fstream>
#include "../../Common/Common.h"
#include "../../Common/Utilities.hpp"
#include "../../Common/Timer.h"
#include "../../Graph/Vertex/CPU/BSTVertex.h"
#include "../../Utility/ExperimentConfiguration.hpp"

/*
 *
 */

#define SerialBA_DEBUG 0
#define SerialBA_SIMULATED 0
#define SerialBA_INFO 0
#define SerialBA_LOG 1
#define SerialBA_OUTPUT_FILE 0
#define SerialBA_OUTPUT_DEGREE_DIST 1
#define SerialBA_RETRY_COUNT 0

class SerialBA
{
private:
    //#Node, #Edge per Node
    NodeType _N;
    unsigned int _M;

    int findNode(const NodeType degree[], const NodeType N, NodeType R);

    // Graph are treated as G[N][d]
    NodeType* initializeNetwork(std::ofstream& log);
    void clearNetwork(NodeType *_G);

    inline void randomSet(NodeType value[], std::set<NodeType>& set, EdgeType Z) __attribute__((always_inline));
    void generateGraphNetworkX(ExperimentConfiguration conf);

public:
    SerialBA(LongLong n, LongLong m = 1);
    virtual ~SerialBA();

    void generateGraph(ExperimentConfiguration conf);
};

#endif /* SERIALBA_ALGO_H_ */
