/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "HashVertex.h"

#include <cstring>

HashVertex::HashVertex(EdgeIndexType n)
{
    _n = n;
    _size = (EdgeIndexType) ((long double) n * OVERHEAD + 1.0);
    _A = new NodeType[_size];
    std::memset(_A, NULL_ID, _size * sizeof(NodeType));
    _count = 0;
}

HashVertex::~HashVertex()
{
    delete[] _A;
    _A = NULL;
}

/* hashing with open addressing by division */
EdgeIndexType HashVertex::insertItem(NodeType item)
{
    if (_count < _n)
    {
        long long probe;
        for (probe = (item * 2654435761) % _size; _A[probe] != NULL_ID; probe = (probe + 1) % _size)
        {
            if (_A[probe] == item)
                return NULL_ID;
        }
        _A[probe] = item;
        _count++;
        return probe;
    }
    else //Full
    {
        return NULL_ID;
    }
}

bool HashVertex::contains(NodeType item)
{
    for (long long probe = (item * 2654435761) % _size; _A[probe] != NULL_ID; probe = (probe + 1) % _size)
    {
        if (_A[probe] == item)
        {
            return true;
        }
    }
    return false;
}

EdgeIndexType HashVertex::indexOf(NodeType item)
{
    for (long long probe = (item * 2654435761) % _size; _A[probe] != NULL_ID; probe = (probe + 1) % _size)
    {
        if (_A[probe] == item)
        {
            return probe;
        }
    }
    return NULL_ID;
}

