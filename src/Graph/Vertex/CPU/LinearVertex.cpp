/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "LinearVertex.h"
#include <iostream>
#include <limits.h>

LinearVertex::LinearVertex(EdgeIndexType d)
{
    this->_root = 0;
    this->_M = d;
    this->_size = d;
    this->_count = 0;
    this->_A = new NodeType[d];
    this->_child = new EdgeIndexType*[2];
    this->_child[LEFT_CHILD] = new EdgeIndexType[d];
    this->_child[RIGHT_CHILD] = new EdgeIndexType[d];

    // -1 is important
    std::fill_n(_A, _M, -1);
    std::fill_n(_child[LEFT_CHILD], _M, NULL_PTR);
    std::fill_n(_child[RIGHT_CHILD], _M, NULL_PTR);
}

LinearVertex::~LinearVertex()
{
    delete[] _child[RIGHT_CHILD];
    delete[] _child[LEFT_CHILD];
    delete[] _child;
    delete[] _A;
}

