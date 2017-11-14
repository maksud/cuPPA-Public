/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "cuRandom.cuh"
#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <curand_kernel.h>
#include <ctime>
#include <cstdio>

__device__ __constant__ unsigned long long CONST_RAND_BIG_MAX = 0x3FFFFFFFFFFFFFFF;
__device__ __constant__ unsigned long long CONST_RAND_BIG_MAX_PLUS_ONE = 0x4000000000000000;
__device__ __constant__ unsigned long long CONST_RAND_MAX_PLUS_ONE = 0x80000000;

__device__ unsigned long long getUniformRandomNN(curandState* state)
{
    return curand_uniform_double(state) * CONST_RAND_BIG_MAX;
}

__device__ unsigned int getUniformRandomN(curandState* state)
{
    return curand_uniform(state) * 4294967295;
}

__device__ double getUniformRandom01(curandState* state)
{
    return curand_uniform_double(state);
}

__device__ curandState initState(unsigned long long seed = 0)
{
    curandState state;
    curand_init((unsigned long long) seed, 0, 0, &state);
    return state;
}

__global__ void setup_curand_kernel_MRG32k3a(curandStateMRG32k3a *state, unsigned long long seed)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    /* Each thread gets same seed, a different sequence number, no offset */
    curand_init(seed, id, 0, &state[id]);
}

__global__ void setup_curand_kernel(curandState *state, unsigned long long seed)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    /* Each thread gets same seed, a different sequence number, no offset */
    curand_init(seed, id, 0, &state[id]);
}

