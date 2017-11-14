/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include <iostream>
#include <numeric>
#include <stdlib.h>

static void CheckCudaErrorAux(const char *, unsigned, const char *, cudaError_t);
#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)
