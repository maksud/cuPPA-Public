/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef _RANDOM64_HPP
#define _RANDOM64_HPP

#include <cstdio>
#include <stdint.h>    // integer types

/* Period parameters */
#define __MT19937_NN 312
#define __MT19937_MM 156
#define __MT19937_MATRIX_AA 0xB5026F5AA96619E9ULL
#define __MT19937_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define __MT19937_LM 0x7FFFFFFFULL /* Least significant 31 bits */

// -- 64 bit
/* The array for the state vector */
static unsigned long long __MT19937_mt[__MT19937_NN];
/* mti==NN+1 means mt[NN] is not initialized */
static int __MT19937_mti = __MT19937_NN + 1;

class Random64
{
public:

    /* initializes mt[NN] with a seed */
    static void Seed(unsigned long long seed)
    {
        __MT19937_mt[0] = seed;
        for (__MT19937_mti = 1; __MT19937_mti < __MT19937_NN; __MT19937_mti++)
        {
            __MT19937_mt[__MT19937_mti] = (6364136223846793005ULL * (__MT19937_mt[__MT19937_mti - 1] ^ (__MT19937_mt[__MT19937_mti - 1] >> 62)) + __MT19937_mti);
        }
    }

    /* generates a random number on [0,1]-real-interval */
    inline static double Uniform01Inclusive(void)
    {
        return (double) (Uniform() >> 11) * (1.0 / 9007199254740991.0);
    }

    /* generates a random number on [0,1)-real-interval */
    inline static double Uniform01(void)
    {
        return (double) (Uniform() >> 11) * (1.0 / 9007199254740992.0);
    }

    /* generates a random number on (0,1)-real-interval */
    inline double Uniform01Exclusive(void)
    {
        return ((double) (Uniform() >> 12) + 0.5) * (1.0 / 4503599627370496.0);
    }

    /* generates a random number on [0, 2^63-1]-interval */
    inline long long Uniform63(void)
    {
        return (long long) (Uniform() >> 1);
    }

    /* generates a random number on [0, 2^64-1]-interval */
    inline static uint64_t Uniform(void)
    {
        int i;
        unsigned long long x;
        static unsigned long long mag01[2] = { 0ULL, __MT19937_MATRIX_AA };

        if (__MT19937_mti >= __MT19937_NN)
        { /* generate NN words at one time */

            /* if Seed64() has not been called, */
            /* a default initial seed is used     */
            if (__MT19937_mti == __MT19937_NN + 1)
                Seed(5489ULL);

            for (i = 0; i < __MT19937_NN - __MT19937_MM; i++)
            {
                x = (__MT19937_mt[i] & __MT19937_UM) | (__MT19937_mt[i + 1] & __MT19937_LM);
                __MT19937_mt[i] = __MT19937_mt[i + __MT19937_MM] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            for (; i < __MT19937_NN - 1; i++)
            {
                x = (__MT19937_mt[i] & __MT19937_UM) | (__MT19937_mt[i + 1] & __MT19937_LM);
                __MT19937_mt[i] = __MT19937_mt[i + (__MT19937_MM - __MT19937_NN)] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            x = (__MT19937_mt[__MT19937_NN - 1] & __MT19937_UM) | (__MT19937_mt[0] & __MT19937_LM);
            __MT19937_mt[__MT19937_NN - 1] = __MT19937_mt[__MT19937_MM - 1] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];

            __MT19937_mti = 0;
        }

        x = __MT19937_mt[__MT19937_mti++];

        x ^= (x >> 29) & 0x5555555555555555ULL;
        x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
        x ^= (x << 37) & 0xFFF7EEE000000000ULL;
        x ^= (x >> 43);

        return x;
    }

    /* generates a random number on [start, end]-interval */
    inline static uint64_t UniformRange(uint64_t start, uint64_t end)
    {
        return (start + (Uniform() % (end - start + 1)));
    }

    /* generates a random number on [0, max-1]-interval */
    inline static uint64_t UniformMax(uint64_t max)
    {
        return (Uniform() % max);
    }

    /* generates a random number on [1, max]-interval */
    inline static uint64_t UniformMax1(LongLong max)
    {
        return ((Uniform() % max) + 1);
    }
};

#endif // _RANDOM64_HPP
