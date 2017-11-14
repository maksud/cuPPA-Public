/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "PAGraph.h"

PAGraph::PAGraph(NodeIndexType n, EdgeIndexType d)
{
    this->_M = d;
    this->_N = n;

    this->_root = new EdgeIndexType[_N];
    this->_size = new EdgeIndexType[_N];
    this->_V = new NodeType*[_N];
    this->_child = new EdgeIndexType**[_N];

    for (NodeType i = 0; i < _N; i++)
    {
        this->_V[i] = new NodeType[_M];
        this->_child[i] = new EdgeIndexType*[2];
        this->_child[i][LEFT_CHILD] = new EdgeIndexType[d];
        this->_child[i][RIGHT_CHILD] = new EdgeIndexType[d];

        fill_n(_V[i], _M, LLONG_MAX);
        fill_n(_child[i][LEFT_CHILD], _M, NULL_PTR);
        fill_n(_child[i][RIGHT_CHILD], _M, NULL_PTR);
    }
}

PAGraph::~PAGraph()
{
}

void PAGraph::Adj(NodeIndexType nodeIndex, EdgeIndexType edgeIndex, NodeType v)
{
    _V[nodeIndex][edgeIndex] = v;
}

NodeType PAGraph::Adj(NodeIndexType nodeIndex, EdgeIndexType edgeIndex)
{
    return _V[nodeIndex][edgeIndex];
}
