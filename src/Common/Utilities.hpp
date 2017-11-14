/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdio>
#include <sys/time.h>
#include "DataTypes.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <unistd.h>
#include <fstream>
#include "Common.h"

#define IS_MAC true

class Utilities
{
public:
    /*-----------------------FindMax------------------------------------------*/
    template<typename T, typename SizeT>
    static T FindMax(const T *d, SizeT size)
    {
        T max = d[0];
        for (SizeT i = 1; i < size; i++)
            if (max < d[i])
                max = d[i];
        return max;
    }

    /*-----------------------FindMin------------------------------------------*/
    template<typename T, typename SizeT>
    static T FindMin(const T *d, SizeT size)
    {
        T min = d[0];
        for (SizeT i = 1; i < size; i++)
            if (min > d[i])
                min = d[i];
        return min;
    }

    static std::string getCurrenttime()
    {
        time_t rawtime;
        struct tm * timeinfo;
        char buffer[80];

        time(&rawtime);
        timeinfo = localtime(&rawtime);

        strftime(buffer, 80, "%d-%m-%Y %I:%M:%S", timeinfo);
        std::string str(buffer);
        return str;
    }

    static long difftime(struct timeval &end_time, struct timeval &start_time)
    {
        struct timeval difference;

        difference.tv_sec = end_time.tv_sec - start_time.tv_sec;
        difference.tv_usec = end_time.tv_usec - start_time.tv_usec;

        /* Using while instead of if below makes the code slightly more robust. */
        while (difference.tv_usec < 0)
        {
            difference.tv_usec += 1000000;
            difference.tv_sec -= 1;
        }
        return (1000000L * difference.tv_sec + difference.tv_usec);
    } /* timeval_diff() */

    static int parseLine(char* line)
    {
        int i = strlen(line);
        while (*line < '0' || *line > '9')
            line++;
        line[i - 3] = '\0';
        i = atoi(line);
        return i;
    }

#if IS_MAC
    static int getUsedVirtualMemory()
    {
        return 0;
    }
#else
    static int getUsedVirtualMemory()
    { //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL)
        {
            if (strncmp(line, "VmSize:", 7) == 0)
            {
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
#endif

#if IS_MAC
    static int getUsedPhysicalMemory()
    {
        return 0;
    }
#else
    static int getUsedPhysicalMemory()
    { //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL)
        {
            if (strncmp(line, "VmRSS:", 6) == 0)
            {
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
#endif

    static std::string getSelfPath()
    {
        char buff[4096];
        ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff) - 1);
        if (len != -1)
        {
            buff[len] = '\0';
            return std::string(buff);
        }
        else
        {
            /* handle error condition */
            return std::string("");
        }
    }

    static void logMemUsage(std::ostream* _log, std::string message)
    {
        (*_log) << std::endl << message << std::endl;
        (*_log) << "-----------------------------------------------" << std::endl;
        (*_log) << message << ":VM: " << getUsedVirtualMemory() << std::endl;
        (*_log) << message << ":PM: " << getUsedPhysicalMemory() << std::endl << std::endl;
    }

    static unsigned int upperPower2(unsigned int x)
    {
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return ++x;
    }

    static unsigned int lowerPower2(unsigned int n)
    {
        int p = 1;
        while (p <= n)
            p <<= 1;
        p >>= 1;
        return p;
    }

};

#endif /* UTILITIES_HPP_ */
