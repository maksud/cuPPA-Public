/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include <stdint.h>

#define NODE32BIT

#ifndef NODE32BIT
typedef int64_t NodeType;
typedef int64_t NodeIndexType;
#else
typedef int32_t NodeType;
typedef int32_t NodeIndexType;
#endif

//For edge
#ifndef EDGE32BIT
typedef int64_t EdgeType;
typedef int16_t EdgeIndexType;
#else
typedef int32_t EdgeType;
typedef int16_t EdgeIndexType;
#endif
//For Processor
typedef int ProcessorType;

//For General
typedef long long LongLong;

//For Node
typedef unsigned int uint;

#ifdef SHORT_WEIGHT
typedef short int WeightType;
typedef short int DegreeType;
#else
typedef unsigned int WeightType;
typedef unsigned int DegreeType;
#endif

typedef double ProbabilityType;

struct Distribution
{
    unsigned int Delta;
    WeightType* Weights;
    NodeType* Frequencies;
};

template<typename Degree, typename Value>
struct DistributionKV
{
    unsigned int Size;
    Degree *Keys;
    Value *Values;
};

#define USE_GRAPH_CONTAINER 0

#endif /* DATATYPES_H_ */
