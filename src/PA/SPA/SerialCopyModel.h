/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef PA_SERIAL_SERIALCOPYMODEL_H_
#define PA_SERIAL_SERIALCOPYMODEL_H_

#include <fstream>
#include <string>
#include <set>
#include "../../Common/Common.h"
#include "../../Common/Utilities.hpp"
#include "../../Common/Timer.h"
#include "../../Graph/Vertex/CPU/BSTVertex.h"
#include "../../Utility/ExperimentConfiguration.hpp"

#define SerialCopyModel_DEBUG 0
#define SerialCopyModel_SIMULATED 0
#define SerialCopyModel_INFO 0
#define SerialCopyModel_LOG 1
#define SerialCopyModel_OUTPUT_FILE 0
#define SerialCopyModel_OUTPUT_DEGREE_DIST 1
#define SerialCopyModel_RETRY_COUNT 1

class SerialCopyModel
{
private:
    //#Node, #Edge per Node
    NodeType _N;
    EdgeIndexType _M;
    double _p;

    int findNode(const NodeType degree[], const NodeType N, NodeType R);

    // Graph are treated as G[N][d]
    NodeType* initializeNetwork(std::ofstream& log, int init_val);
    void clearNetwork(NodeType *_G);
    void generateGraphCopyModelBST(ExperimentConfiguration conf);

public:
    SerialCopyModel(NodeType n, EdgeIndexType m, double p);
    virtual ~SerialCopyModel();

    void generateGraph(ExperimentConfiguration conf);

};

#endif /* PA_SERIAL_SERIALCOPYMODEL_H_ */
