/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef _RANDOM32_HPP
#define _RANDOM32_HPP

#include <cstdio>
#include <stdint.h>    // integer types

/* Period parameters */
#define __MT19937_N 624
#define __MT19937_M 397
#define __MT19937_MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define __MT19937_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define __MT19937_LOWER_MASK 0x7fffffffUL /* least significant r bits */

// -- 32 bit
/* The array for the state vector */
static unsigned long mt[__MT19937_N];
/* mti==N+1 means mt[N] is not initialized */
static int mti = __MT19937_N + 1;

class Random32
{
public:

    /* initializes mt[N] with a seed */
    static void Seed(unsigned long seed)
    {
        mt[0] = seed & 0xffffffffUL;
        for (mti = 1; mti < __MT19937_N; mti++)
        {
            mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
            /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
            /* In the previous versions, MSBs of the seed affect   */
            /* only MSBs of the array mt[].                        */
            /* 2002/01/09 modified by Makoto Matsumoto             */
            mt[mti] &= 0xffffffffUL;
            /* for >32 bit machines */
        }
    }

    /* generates a random number on [0,1]-real-interval */
    inline static double Uniform01Inclusive(void)
    {
        return Uniform() * (1.0 / 4294967295.0);
    }

    /* generates a random number on [0,1)-real-interval */
    inline static double Uniform01(void)
    {
        return Uniform() * (1.0 / 4294967296.0);
    }

    /* generates a random number on (0,1)-real-interval */
    inline double Uniform01Exclusive(void)
    {
        return (((double) Uniform()) + 0.5) * (1.0 / 4294967296.0);
    }

    /* generates a random number on [0, 2^63-1]-interval */
    inline long long Uniform63(void)
    {
        return (long) (Uniform() >> 1);
    }

    /* generates a random number on [0,0xffffffff]-interval */
    inline static uint32_t Uniform(void)
    {
        unsigned long y;
        static unsigned long mag01[2] = { 0x0UL, __MT19937_MATRIX_A };
        /* mag01[x] = x * MATRIX_A  for x=0,1 */

        if (mti >= __MT19937_N)
        { /* generate N words at one time */
            int kk;

            /* if Seed32() has not been called, */
            /* a default initial seed is used */
            if (mti == __MT19937_N + 1)
                Seed(5489UL);

            for (kk = 0; kk < __MT19937_N - __MT19937_M; kk++)
            {
                y = (mt[kk] & __MT19937_UPPER_MASK) | (mt[kk + 1] & __MT19937_LOWER_MASK);
                mt[kk] = mt[kk + __MT19937_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            for (; kk < __MT19937_N - 1; kk++)
            {
                y = (mt[kk] & __MT19937_UPPER_MASK) | (mt[kk + 1] & __MT19937_LOWER_MASK);
                mt[kk] = mt[kk + (__MT19937_M - __MT19937_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            y = (mt[__MT19937_N - 1] & __MT19937_UPPER_MASK) | (mt[0] & __MT19937_LOWER_MASK);
            mt[__MT19937_N - 1] = mt[__MT19937_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

            mti = 0;
        }

        y = mt[mti++];

        /* Tempering */
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        return y;
    }

    /* generates a random number on [start, end]-interval */
    inline static uint32_t UniformRange(uint32_t start, uint32_t end)
    {
        return (start + (Uniform() % (end - start + 1)));
    }

    /* generates a random number on [0, max-1]-interval */
    inline static uint32_t UniformMax(uint32_t max)
    {
        return (Uniform() % max);
    }

    /* generates a random number on [1, max]-interval */
    inline static uint32_t UniformMax1(uint32_t max)
    {
        return ((Uniform() % max) + 1);
    }
};

#endif // _RANDOM32_HPP
