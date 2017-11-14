/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef COMMON_HELPER_CUH_
#define COMMON_HELPER_CUH_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "Common.h"

using namespace std;

#define NODEBUG 1
#define USECUDA 1

#define THREADS_PER_BLOCK 32

/*----------------------------------------------------------------------------*/
extern int _dbglvl;
//extern __device__ int _cudbglvl;

#if NODEBUG
#define _DBG( _lvl, _msg, _el ) ((void)0)
#else
#define _DBG( _lvl, _msg, _el ) if( (_lvl) > _dbglvl ) /*No-op*/; else{\
    if( (_lvl) > 4 ) {\
        cout << __FILE__ << ":" << __LINE__ << ": ";\
    } \
    cout << _msg; \
    if( _el ) cout << endl; \
  }
#endif /*NODEBUG*/
#define DBGL( _lvl, _msg ) _DBG( _lvl, _msg, true )
#define DBG( _lvl, _msg ) _DBG( _lvl, _msg, false )

#define eprintf(...) fprintf (stderr, __VA_ARGS__)

#if NODEBUG
#define cuDBG( _culvl, ... ) ((void)0)
#else
#define cuDBG( _lvl, ... ) if( (_lvl) > _cudbglvl ) /*No-op*/; else{\
    if( (_lvl) > 4 ) {\
        printf("%s\t%d: ", __FILE__, __LINE__);\
    } \
    printf(__VA_ARGS__); \
  }
#endif /*NODEBUG*/

#if NODEBUG
#define cuDBG0( _culvl, ... ) ((void)0)
#else
#define cuDBG0( _lvl, ... ) if( (_lvl) > _cudbglvl ) /*No-op*/; else{\
    if(blockIdx.x * blockDim.x + threadIdx.x == 0){\
        if( (_lvl) > 4 ) {\
            printf("%s\t%d: ", __FILE__, __LINE__);\
        }\
        printf(__VA_ARGS__);\
    }\
  }
#endif /*NODEBUG*/

/*----------------------------------------------------------------------------*/
extern int _enslvl;
#if NODEBUG
#define ENSURE( _lvl, _cond, _msg ) ((void)0)
#define FAIL( _msg ) ((void)0)
#else
#define ENSURE( _lvl, _cond, _msg ) if( (_lvl) > _enslvl ) /*No-op*/; else{\
    if( !( _cond ) ) {\
        cerr << __FILE__ << ":" << __LINE__ << " failed to ensure" << endl;\
        cerr << "\"" << #_cond << "\"" << endl; \
        cerr << "\"" << _msg << "\"" << endl; \
        cerr << "Bailing out..." << endl; \
        exit(1); \
    } \
  }
#define FAIL( _msg ) ENSURE( 0, 0, _msg )
#endif /*NODEBUG*/

/*----------------------------------------------------------------------------*/
#define SQR(_) ((_)*(_))
#define MAX(_a,_b) ((_a)>(_b)?(_a):(_b))
#define SWAP(_t,_a,_b) do{_t _tv = _a; _a = _b; _b = _tv;}while(0)

/*----------------------------------------------------------------------------*/
#include <sys/time.h>
typedef struct timeval TIMER_TYPE;

#define TIMER_NOW(_t) gettimeofday(&_t,NULL)
#define TIMER_US(_t) ((double)(_t).tv_sec*1e6 + (_t).tv_usec)
#define TIMER_DIFF(_t2, _t1) (TIMER_US(_t2)-TIMER_US(_t1))

double WCSEC(TIMER_TYPE t1);

/*----------------------------------------------------------------------------*/
#define ENV_PBOOL(_v,_s) cout<<_s<<"="<<((_v)?"TRUE":"FALSE")<<endl
#define ENV_PRINT(_v,_s) cout<<_s<<"="<<_v<<endl
#define ENV_BOOL(_v,_s,_i) _v = (!getenv(_s) ? (_i) : \
                                 (!strcmp(getenv(_s),"TRUE") || \
                                  !strcmp(getenv(_s),"true") || \
                                  !strcmp(getenv(_s),"1"))); ENV_PBOOL(_v,_s)
#define ENV_DBOOL(_v,_s,_i) bool ENV_BOOL(_v,_s,_i)
#define ENV_INT(_v,_s,_i) _v = (!getenv(_s) ? (_i) : atoi(getenv(_s))); ENV_PRINT(_v,_s)
#define ENV_DINT(_v,_s,_i) int ENV_INT(_v,_s,_i)
#define ENV_DBL(_v,_s,_i) _v = (!getenv(_s) ? (_i) : atof(getenv(_s))); ENV_PRINT(_v,_s)
#define ENV_DDBL(_v,_s,_i) double ENV_DBL(_v,_s,_i)
#define ENV_STR(_v,_s,_i) _v = (!getenv(_s) ? (_i) : getenv(_s)); ENV_PRINT(_v,_s)
#define ENV_DSTR(_v,_s,_i) const char *ENV_STR(_v,_s,_i)

/*----------------------------------------------------------------------------*/
#define EPS 1e-3lf
#define ABS(_) ((_)<0?(-(_)):(_))
#define EQUAL( a, b ) (ABS((a)-(b)) <= EPS)

/*----------------------------------------------------------------------------*/
#define HOSTCalloc(_A, _B, _C) do{ \
        (*(_A)) = (_C *)calloc((_B), sizeof(_C)); \
        ENSURE( 0, (*(_A)), "" ); \
    }while (0)
#define HOSTCalloc2D(_A, _B, _C, _D) do{ \
        _D **_ar = (_D **)calloc((_B), sizeof(_D *)); \
        ENSURE( 0, _ar, "" ); \
        _D *_el = (_D *)calloc((_B)*(_C), sizeof(_D)); \
        ENSURE( 0, _el, "" ); \
        for( int i = 0; i < (_B); i++ ) \
        { \
            _ar[i] = _el + i*(_C); \
        } \
        (*(_A)) = _ar; \
    }while (0)
#define HOSTFree(_A) do{ \
        if(!(_A)) break; \
        free(_A); \
        (_A) = 0; \
    }while (0)
#define HOSTFree2D(_A) do{ \
        if(!(_A)) break; \
        free( (_A)[0] ); \
        free(_A); \
        (_A) = 0; \
    }while (0)

/*----------------------------------------------------------------------------*/
#define CUDAIAMT0() (threadIdx.x == 0 && threadIdx.y == 0)
#define CUDACALL(_C) do{ \
        cudaError_t _cudaerr = _C; \
        if(_cudaerr != cudaSuccess) \
        { \
            cerr<<__FILE__<<":"<<__LINE__<<": "; \
            cerr<<"CUDACALL failed \""<<(#_C)<<"\": returned " \
                <<_cudaerr<<endl; \
            exit(1); \
        } \
    }while (0)
#define CUDAMemcpyH2D(_A, _B, _C, _D) do{ \
        CUDACALL(cudaMemcpy((_A), (_B), (_C)*sizeof(_D), cudaMemcpyHostToDevice)); \
    }while (0)
#define CUDAMemcpyD2H(_A, _B, _C, _D) do{ \
        CUDACALL(cudaMemcpy((_A), (_B), (_C)*sizeof(_D), cudaMemcpyDeviceToHost)); \
    }while (0)
#define CUDAMalloc(_A, _B) do{ \
        CUDACALL(cudaMalloc((_A), (_B))); \
    }while (0)
#define CUDAMemset(_A, _v, _B, _C) do{ \
        CUDACALL(cudaMemset((_A), _v, (_B)*sizeof(_C))); \
    }while (0)
#define CUDACalloc(_A, _B, _C) do{ \
        CUDAMalloc((_A), (_B)*sizeof(_C)); \
        CUDAMemset((*(_A)), 0, (_B), _C); \
    }while (0)
#define CUDAAllocSet(_A, _B, _C, _D) do{ \
        CUDAMalloc((_A), (_C)*sizeof(_D)); \
        CUDAMemcpyH2D((*(_A)), (_B), (_C), _D); \
    }while (0)
#define CUDAFree(_A) cudaFree(_A)
#define CUDAAtomicAdd(_A,_B) atomicAdd((_A),(_B))
#define CUDAAtomicSub(_A,_B) atomicSub((_A),(_B))

#if USECUDA
#define HOSTPinnedAlloc(_A, _B, _C) do{ \
        CUDACALL(cudaMallocHost((_A), (_B)*sizeof(_C))); \
    }while (0)
#define HOSTPinnedFree(_A) CUDACALL(cudaFreeHost(_A))
#else/*USECUDA*/
#define HOSTPinnedAlloc(_A, _B, _C) HOSTCalloc(_A,_B,_C)
#define HOSTPinnedFree(_A) HOSTFree(_A)
#endif/*USECUDA*/

int getSMCount();

#endif /* COMMON_HELPER_CUH_ */
