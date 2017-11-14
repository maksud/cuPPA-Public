/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef CUPPA_CUH_
#define CUPPA_CUH_

#define cuPPA_DEBUG 0
#define DEBUG_cuPPA_RRP 0

#define USE_DEFAULT_CURAND 1 //Default CURAND performs well

#define SHOW_INTERNAL_TIMINGS 1
#define SHOW_QUEUE_OUTPUT 0
#define WAITING_QUEUE_IN_GLOBAL_MEMORY 0

#include "../../Common/DataTypes.h"
#include "../../Common/CommonHelper.cuh"
#include "../../Utility/ExperimentConfiguration.hpp"
#include <curand_kernel.h>
#include "../../Common/cuRandom.cuh"
#include "../../Common/Timer.h"
#include "../../Common/Utilities.hpp"
#include "../../Graph/Vertex/CUDA/cuBSTVertex.cuh"
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>

#if USE_DEFAULT_CURAND
__device__ bool re_execute_copy_model(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandState *state);
#else
__device__ bool re_execute_copy_model(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandStateMRG32k3a *state);
#endif

__global__ void initialize_vertex_structures(NodeType nVertices, EdgeIndexType M, NodeType *A, EdgeIndexType *children, cuBSTVertex *vertices);
__global__ void initial_graph(EdgeIndexType M, cuBSTVertex *vertices);

#if USE_DEFAULT_CURAND
__global__ void generate_graph_race(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandState *state);
__global__ void generate_graph_mutex(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandState *state);
__global__ void generate_graph_avoid(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandState *state);
#else
__global__ void generate_graph_race(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state);
__global__ void generate_graph_mutex(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state);
__global__ void generate_graph_avoid(NodeType nVertices, EdgeIndexType M, double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state);
#endif

//
__global__ void find_max_queue_size(int X[]);

//Blocked RRP Specific Implementation
__global__ void validate1_graph_rrp_by_batch(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F);
__global__ void validate2_graph_rrp_by_batch(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F);
__global__ void initialize_vertex_structures_rrp(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F);
__global__ void initialize_vertex_structures_rrp_fast(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F);
__global__ void initialize_vertex_structures_rrp_2(NodeType u, NodeType B, EdgeIndexType M, NodeType *G, EdgeIndexType *children, cuBSTVertex *F);
__global__ void initial_graph_rrp(EdgeIndexType M, NodeType *G);

#if USE_DEFAULT_CURAND
__device__ bool re_execute_copy_model_rrp(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandState *state);
#else
__device__ bool re_execute_copy_model_rrp(NodeType *G, cuBSTVertex *V, EdgeIndexType M, double p, NodeType u, EdgeIndexType i, curandStateMRG32k3a *state);
#endif

#if USE_DEFAULT_CURAND
__global__ void generate_graph_rrp_by_batch(const NodeType u, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandState *state, int *MaxQSize);
#else
__global__ void generate_graph_rrp_by_batch(const NodeType u, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state, int *MaxQSize);
#endif

#if USE_DEFAULT_CURAND
__global__ void generate_graph_rrp_by_batch_execute_copy_model(const NodeType u, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandState *state, EdgeType* Q, int *MaxQSize);
__global__ void generate_graph_rrp_by_batch_resolve_queue(const NodeType u, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandState *state, EdgeType* Q, int *MaxQSize);
#else
__global__ void generate_graph_rrp_by_batch_execute_copy_model(const NodeType u, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state, EdgeType* Q, int *MaxQSize);
__global__ void generate_graph_rrp_by_batch_resolve_queue(const NodeType u, const NodeType B, const EdgeIndexType M, const double p, cuBSTVertex *F, NodeType *G, curandStateMRG32k3a *state, EdgeType* Q, int *MaxQSize);
#endif

class cuPPA
{
protected:
    //#Node, #Edge per Node
    NodeType _N;
    EdgeIndexType _M;
    double _p;

    NodeType *_dGraph;

public:
    cuPPA(NodeType nNodes, EdgeIndexType mEdgesPerNode, double p);
    virtual ~cuPPA();

    void initialize();

    void generateGraphOld(ExperimentConfiguration conf);
    void generateGraph(ExperimentConfiguration conf);
    void generateGraphAdaptive(ExperimentConfiguration conf);

    void printDegreeDistribution(ExperimentConfiguration conf);
};

#endif /* CUPPA_CUH_ */
