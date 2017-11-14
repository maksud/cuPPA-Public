/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef HASHVERTEX_H_
#define HASHVERTEX_H_

#include "Vertex.h"

#define MIN_ID (0)
#define MAX_ID (999999999)

/* overhead parameter that determines both space and search costs */
/* must be strictly greater than 1 */
#define OVERHEAD (1)
#define NULL_ID (-1)

class HashVertex: public Vertex
{
    EdgeIndexType _n;
    EdgeIndexType _size;
    NodeType* _A;
    EdgeIndexType _count;

public:
    HashVertex(EdgeIndexType n);
    virtual ~HashVertex();

    EdgeIndexType size()
    {
        return _size;
    }

    EdgeIndexType count()
    {
        return _count;
    }

    EdgeIndexType insertItem(NodeType v);
    bool contains(NodeType v);
    EdgeIndexType indexOf(NodeType v);
};

#endif /* HASHVERTEX_H_ */
