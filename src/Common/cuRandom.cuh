/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef CUDARANDOM_CUH_
#define CUDARANDOM_CUH_

#include <iostream>
#include "cuda_runtime.h"
#include <device_launch_parameters.h>
#include <curand_kernel.h>
#include <ctime>
#include <cstdio>

__device__ unsigned long long getUniformRandomNN(curandState* state);
__device__ unsigned int getUniformRandomN(curandState* state);
__device__ double getUniformRandom01(curandState* state);
__device__ curandState initState(unsigned long long seed);

__global__ void setup_curand_kernel_MRG32k3a(curandStateMRG32k3a *state, unsigned long long seed);
__global__ void setup_curand_kernel(curandState *state, unsigned long long seed);

#endif /* CUDARANDOM_CUH_ */

