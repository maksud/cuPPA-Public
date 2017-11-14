/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#include "CommonHelper.cuh"

int _dbglvl = 5;
int _enslvl = 10;
__device__ int _cudbglvl = 5;

double WCSEC(TIMER_TYPE t1)
{
    TIMER_TYPE t2;
    TIMER_NOW(t2);
    return TIMER_DIFF(t2,t1) / 1e6;
}

int getSMCount(void)
{
    static int nprocs = -1;
    if (nprocs < 0)
    {
        cudaDeviceProp cdprop;
        CUDACALL(cudaGetDeviceProperties(&cdprop, 0));
        nprocs = cdprop.multiProcessorCount;
        ENSURE(0, nprocs > 0, "");
    }
    return nprocs;
}
