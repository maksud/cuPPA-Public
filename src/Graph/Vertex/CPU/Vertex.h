/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef VERTEX_H_
#define VERTEX_H_

#include "../../../Common/DataTypes.h"
#include "../../../Common/Common.h"

class Vertex
{
public:
    virtual ~Vertex()
    {
    }
    virtual EdgeIndexType size() = 0;

    virtual EdgeIndexType count() = 0;

    virtual EdgeIndexType insertItem(NodeType v) = 0;

    virtual bool contains(NodeType v) = 0;

    virtual EdgeIndexType indexOf(NodeType v) = 0;
};

#endif /* VERTEX_H_ */
