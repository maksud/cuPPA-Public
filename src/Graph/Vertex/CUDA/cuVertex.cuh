/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef CUVERTEX_H_
#define CUVERTEX_H_

#include "../../../Common/DataTypes.h"
#include "../../../Common/Common.h"

class cuVertex
{
public:

    CUDA_CALLABLE_MEMBER EdgeIndexType size();

    CUDA_CALLABLE_MEMBER EdgeIndexType count();

    CUDA_CALLABLE_MEMBER bool insertItem(EdgeIndexType i, NodeType v);

    CUDA_CALLABLE_MEMBER bool contains(NodeType v);

    CUDA_CALLABLE_MEMBER EdgeIndexType indexOf(NodeType v);
};

#endif /* VERTEX_H_ */
