/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <stdint.h>

using namespace std;

#define RAND_BIG_MAX          0x3FFFFFFFFFFFFFFF
#define DOUBLE_MAX            0x3FFFFFFFFFFFFFFF
#define RAND_BIG_MAX_PLUS_ONE 0x4000000000000000
#define RAND_MAX_PLUS_ONE     0x80000000

#include "Random32.hpp"
#include "Random64.hpp"

#define RAND_64_BIG_MAX          0xFFFFFFFFFFFFFFFF
#define RAND_32_BIG_MAX          0xFFFFFFFF

#define USE_BUILT_IN_RANDOM 0

#if USE_BUILT_IN_RANDOM

class Random
{
public:

    static void Seed(unsigned long long seed)
    {
        srand(seed);
    }

    static void print()
    {
        cout << "RAND_MAX: " << RAND_MAX << endl;
        unsigned long long rand_max = RAND_MAX;
        cout << "rand_max 1: " << rand_max << endl;
        rand_max = rand_max << 31;
        cout << "rand_max 2: " << rand_max << endl;
        rand_max = rand_max | (unsigned long long) RAND_MAX;
        cout << "rand_max 3: " << rand_max << endl;

        double r = (double) (((LongLong) RAND_MAX << 31) | RAND_MAX);
        cout << "r\t" << r << endl;
        cout << "RAND_BIG_MAX: " << RAND_BIG_MAX << endl;
        cout << "RAND_BIG_MAX_PLUS_ONE: " << RAND_BIG_MAX_PLUS_ONE << endl;

        cout << sizeof(double) << " " << sizeof(float) << " " << sizeof(long long) << endl;
        cout << RAND_MAX << endl;

        long r2 = ((LongLong) RAND_MAX << 31) | RAND_MAX;
        cout << "r2: " << r2 << endl;

        double a = 1.0;
        a = a / RAND_BIG_MAX;
        cout << "a\t" << a << endl;
    }

    /* generates a random number on [0,1]-real-interval */
    inline static double Uniform01Inclusive()
    {
        return (Uniform() / RAND_BIG_MAX);
    }

    /* generates a random number on [0,1)-real-interval */
    inline static double Uniform01(void)
    {
        return (Uniform() / RAND_BIG_MAX_PLUS_ONE);
    }

    //  Return uniformly distributed double in range [0.0, 1.0).
    //inline static double Uniform01_32()
    //{
    //	return ((double) rand() / RAND_MAX_PLUS_ONE);
    //}

    //  Return uniformly distributed double in range [0.0, 1.0].
    //inline static double Uniform01NZ_32()
    //{
    //	return ((double) rand() / RAND_MAX);
    //}

    /* generates a random number on [0, 2^62-1]-interval */
    inline static uint64_t Uniform()
    {
        return (((uint64_t) rand() << 31) | rand());
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
#else
class Random: public Random64
{

};
#endif
#endif // _RANDOM_HPP
