/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef CONFIGURATION_HPP_
#define CONFIGURATION_HPP_

#include "ini/minIni.h"
#include "../Common/DataTypes.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <vector>
using namespace std;

#define sizearray(a)  (sizeof(a) / sizeof((a)[0]))

class InputParser
{
public:
    InputParser(int &argc, char **argv)
    {
        for (int i = 1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }

    std::string getOption(const std::string &option, const std::string defaultValue = "")
    {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end())
        {
            return *itr;
        }
        return defaultValue;
    }

    bool hasOption(const std::string &option)
    {
        return std::find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
    }
private:
    std::vector<std::string> tokens;
};

struct ExperimentConfiguration
{
    string _algorithm;

    NodeType _n1;
    WeightType _x1;
    WeightType _d1min;
    WeightType _d1max;
    WeightType _d1avg;
    string _weight1;
    string _distfile1;
    string _avgcc1;
    string _probabilityMatrixFile1;

    NodeType _n2;
    WeightType _x2;
    WeightType _d2min;
    WeightType _d2max;
    WeightType _d2avg;
    string _weight2;
    string _distfile2;
    string _avgcc2;
    string _probabilityMatrixFile2;

    //MPI
    long _msgBufferSize;
    long _mpiBufferSize;

    //CUDA
    int _blockSize;
    int _threadsPerBlock;
    int _warpSize;

    //LOAD BALANCING
    long _stages;
    double _lbFactor;
    double _alpha;
    double _gamma;
    double _probability;
    string _loadBalancing;
    //
    string _outputDir;
    int _write;
    //
    bool _FLAG_SHOW_DEGREE_DISTRIBUTION;
    int _N_STAGE_TIMINGS;

    void print()
    {
        printf("algorithm: \t%s\n", _algorithm.c_str());

        printf("n1: \t%ld\n", _n1);
        printf("x1: \t%d\n", _x1);
        printf("d1min: \t%d\n", _d1min);
        printf("d1max: \t%d\n", _d1max);
        printf("d1avg: \t%d\n", _d1avg);
        printf("algorithm: \t%s\n", _algorithm.c_str());
        printf("algorithm: \t%s\n", _algorithm.c_str());
        printf("algorithm: \t%s\n", _algorithm.c_str());
        printf("algorithm: \t%s\n", _algorithm.c_str());

    }

    void load(const char* inifile)
    {
        char str[512];

        ini_gets("config", "algorithm", "", str, sizearray(str), inifile);
        _algorithm = string(str);

        {
            ini_gets("config", "n1", "1000", str, sizearray(str), inifile);
            _n1 = atoll(str);
            ini_gets("config", "x1", "2", str, sizearray(str), inifile);
            _x1 = atoi(str);
            ini_gets("config", "d1min", "2", str, sizearray(str), inifile);
            _d1min = atoi(str);
            ini_gets("config", "d1max", "10", str, sizearray(str), inifile);
            _d1max = atoi(str);
            ini_gets("config", "d1avg", "5", str, sizearray(str), inifile);
            _d1avg = atoi(str);
            ini_gets("config", "weight1", "none", str, sizearray(str), inifile);
            _weight1 = string(str);
            ini_gets("config", "distfile1", "", str, sizearray(str), inifile);
            _distfile1 = string(str);
            ini_gets("config", "avgccfile1", "", str, sizearray(str), inifile);
            _avgcc1 = string(str);
            ini_gets("config", "pmatfile1", "", str, sizearray(str), inifile);
            _probabilityMatrixFile1 = string(str);
            ini_gets("config", "alpha", "0.5", str, sizearray(str), inifile);
            _alpha = atof(str);
            ini_gets("config", "gamma", "2.7", str, sizearray(str), inifile);
            _gamma = atof(str);
        }

        {
            ini_gets("config", "n2", "0", str, sizearray(str), inifile);
            _n2 = atoll(str);
            ini_gets("config", "x2", "2", str, sizearray(str), inifile);
            _x2 = atoi(str);
            ini_gets("config", "d2min", "2", str, sizearray(str), inifile);
            _d2min = atoi(str);
            ini_gets("config", "d2max", "10", str, sizearray(str), inifile);
            _d2max = atoi(str);
            ini_gets("config", "d2avg", "5", str, sizearray(str), inifile);
            _d2avg = atoi(str);
            ini_gets("config", "weight2", "none", str, sizearray(str), inifile);
            _weight2 = string(str);
            ini_gets("config", "distfile2", "", str, sizearray(str), inifile);
            _distfile2 = string(str);
            ini_gets("config", "avgccfile2", "", str, sizearray(str), inifile);
            _avgcc2 = string(str);
            ini_gets("config", "pmatfile2", "", str, sizearray(str), inifile);
            _probabilityMatrixFile2 = string(str);
        }

        {
            ini_gets("config", "outputdir", "./out/", str, sizearray(str), inifile);
            _outputDir = string(str);
            ini_gets("config", "msg_buffer", "1024", str, sizearray(str), inifile);
            _msgBufferSize = atoi(str);
            ini_gets("config", "mpi_buffer", "1024", str, sizearray(str), inifile);
            _mpiBufferSize = atoi(str);
        }

        {
            ini_gets("config", "lbfactor", "1.0", str, sizearray(str), inifile);
            _lbFactor = atof(str);
            ini_gets("config", "probability", "0.1", str, sizearray(str), inifile);
            _probability = atof(str);

            ini_gets("config", "stages", "1000", str, sizearray(str), inifile);
            _stages = atof(str);
            ini_gets("config", "loadbalancing", "none", str, sizearray(str), inifile);
            _loadBalancing = string(str);
            ini_gets("config", "write", "0", str, sizearray(str), inifile);
            _write = atoi(str);
        }
    }
};

//

enum LoadType
{
    /*Linear Weights*/
    LINEAR_LOAD,
    /*Power-Law Weights*/
    POWER_LAW_LOAD,
    /*Log-Normal Weights*/
    LOG_NORMAL_LOAD,
    /*Constant Weights*/
    CONSTANT_LOAD,
    /*Degree Distribution from File*/
    DEGREE_DIST_FILE,
    /* Not Used */
    NONE
};

enum GenerationAlgorithm
{
    //Barabasi Albert
    METHOD_SERIAL_BA_NETWORKX,
    METHOD_SERIAL_COPY_MODEL,
    // m = 1
    METHOD_PARALLEL_BA_UNIFORM,
    METHOD_PARALLEL_BA_GEOM,
    // General
//	METHOD_PARALLEL_BA_GENERIC_UNIFORM,
//	METHOD_PARALLEL_BA_LINEAR,
//	METHOD_PARALLEL_BA_GENERIC_GEOM,
//	METHOD_PARALLEL_BA_GENERIC_ROUND_ROBIN_OLD,
//	METHOD_PARALLEL_BA_GENERIC_ROUND_ROBIN,
//	METHOD_PARALLEL_BA_GENERIC_ROUND_ROBIN_NEW,
//	METHOD_PARALLEL_BA_GENERIC_THREADED,
    METHOD_PARALLEL_PREFERENTIAL_ATTACHMENT,
    /**/
    METHOD_SERIAL_SMALL_WORLD_NX,
    /**/
    METHOD_SERIAL_SMALL_WORLD,
    /**/
    METHOD_PARALLEL_SMALL_WORLD,

    /**/
    METHOD_9,
    TEST_1,

    /**/
    METHOD_SCL,
    METHOD_FCL,
    METHOD_NCL,
    METHOD_SCL_NDSSL,
    METHOD_SCL_BP,

    /*Parallel Chung-Lu, General Version*/
    METHOD_PCL,

    /*Parallel Chung-Lu, General Version, Limited Memory*/
    METHOD_PCL_LIM_MEM,

    /*Parallel Chung-Lu, NDSSL*/
    METHOD_PCL_NDSSL,

    /*Parallel ER*/
    METHOD_PER,
    METHOD_SER,

    /* BTER */
    METHOD_SERIAL_BTER,
    METHOD_PARALLEL_BTER,

    /* SBM */
    METHOD_SERIAL_SBM,
    METHOD_PARALLEL_SBM,

    /* cuPPA */
    METHOD_CUPPA,
    METHOD_CUPPA_ADAPTIVE
};

enum LoadBalancingScheme
{
    LB_UNIFORM_COST_DIVISIBLE,
    //
    LB_UNIFORM_COST,
    //
    LB_NAIVE,
    LB_LCP,
    LB_BLOCK_CYCLIC_NAIVE,
    LB_BLOCK_CYCLIC_LCP,
    LB_GEOM,
    //
    LB_UNIFORM_EDGE,
    LB_UNIFORM_EDGE_DIVISIBLE,
    //
    LB_RRP,
    LB_RRP_MEM,
    //
    LB_NO_BALANCING,
    //

    LB_MEM_STAGED,
    LB_MEM_STAGED_2P,
    LB_MEM_NON_STAGED,
    LB_MEM_NON_STAGED_2P,
    LB_MEM_BUFFERED
};

class Methods
{
private:
    typedef map<string, GenerationAlgorithm> EnumMapAlgorithm;
    EnumMapAlgorithm enumMapAlgorithm;

    typedef map<string, LoadType> EnumMapLoad;
    EnumMapLoad enumMapLoad;

    typedef map<string, LoadBalancingScheme> EnumMapLoadBalancingScheme;
    EnumMapLoadBalancingScheme enumMapLoadBalancingScheme;
public:

    Methods()
    {
        // Preferential Attachment
        enumMapAlgorithm["sba"] = METHOD_SERIAL_BA_NETWORKX;
        enumMapAlgorithm["scm"] = METHOD_SERIAL_COPY_MODEL;
        enumMapAlgorithm["ppa"] = METHOD_PARALLEL_PREFERENTIAL_ATTACHMENT;

        // Small-World
        enumMapAlgorithm["ssw-nx"] = METHOD_SERIAL_SMALL_WORLD_NX;
        enumMapAlgorithm["ssw"] = METHOD_SERIAL_SMALL_WORLD;
        enumMapAlgorithm["psw"] = METHOD_PARALLEL_SMALL_WORLD;

        enumMapAlgorithm["ncl"] = METHOD_NCL;
        // Miller and Hagberg
        enumMapAlgorithm["scl"] = METHOD_SCL;

        //Pfeiffer III et al.
        enumMapAlgorithm["fcl"] = METHOD_FCL;

        // Maksud and Maleq
        enumMapAlgorithm["scl-ndssl"] = METHOD_SCL_NDSSL;
        enumMapAlgorithm["scl-bp"] = METHOD_SCL_BP;

        // Parallelization of Miller and Hagberg
        enumMapAlgorithm["pcl"] = METHOD_PCL;
        enumMapAlgorithm["pcl-mem"] = METHOD_PCL_LIM_MEM;

        // Parallelization of NDSSL Version
        enumMapAlgorithm["pcl-ndssl"] = METHOD_PCL_NDSSL;

        // ER
        enumMapAlgorithm["ser"] = METHOD_SER;
        enumMapAlgorithm["per"] = METHOD_PER;

        // BTER
        enumMapAlgorithm["sbter"] = METHOD_SERIAL_BTER;
        enumMapAlgorithm["pbter"] = METHOD_PARALLEL_BTER;

        // Stochastic Blockmodels
        enumMapAlgorithm["ssbm"] = METHOD_SERIAL_SBM;
        enumMapAlgorithm["psbm"] = METHOD_PARALLEL_SBM;

        // cuPPA: GPU Preferential-Attachment
        enumMapAlgorithm["cuppa"] = METHOD_CUPPA;
        enumMapAlgorithm["cuppa-adaptive"] = METHOD_CUPPA_ADAPTIVE;

        /**
         * Weight Type
         */
        enumMapLoad["lin"] = LINEAR_LOAD;
        enumMapLoad["pow"] = POWER_LAW_LOAD;
        enumMapLoad["const"] = CONSTANT_LOAD;
        enumMapLoad["none"] = NONE;
        enumMapLoad["distfile"] = DEGREE_DIST_FILE;

        /**
         * Load Balancing Schemes
         */
        enumMapLoadBalancingScheme["unp"] = LB_NAIVE;
        enumMapLoadBalancingScheme["uep"] = LB_UNIFORM_EDGE;
        enumMapLoadBalancingScheme["uep-div"] = LB_UNIFORM_EDGE_DIVISIBLE;
        enumMapLoadBalancingScheme["ucp"] = LB_UNIFORM_COST;
        enumMapLoadBalancingScheme["ucp-div"] = LB_UNIFORM_COST_DIVISIBLE;
        enumMapLoadBalancingScheme["rrp"] = LB_RRP;
        enumMapLoadBalancingScheme["rrp-mem"] = LB_RRP_MEM;
        enumMapLoadBalancingScheme["lcp"] = LB_LCP;
        enumMapLoadBalancingScheme["bcp"] = LB_BLOCK_CYCLIC_NAIVE;
        enumMapLoadBalancingScheme["bcp-lcp"] = LB_BLOCK_CYCLIC_LCP;
        enumMapLoadBalancingScheme["geom"] = LB_GEOM;
        enumMapLoadBalancingScheme["none"] = LB_NO_BALANCING;

        enumMapLoadBalancingScheme["mem-stage"] = LB_MEM_STAGED;
        enumMapLoadBalancingScheme["mem-stage-2p"] = LB_MEM_STAGED_2P;
        enumMapLoadBalancingScheme["mem"] = LB_MEM_NON_STAGED;
        enumMapLoadBalancingScheme["mem-2p"] = LB_MEM_NON_STAGED_2P;
        enumMapLoadBalancingScheme["mem-buf"] = LB_MEM_BUFFERED;
    }

    GenerationAlgorithm ParseAlgorithm(const string& value)
    {
        EnumMapAlgorithm::const_iterator iValue = enumMapAlgorithm.find(value);
        if (iValue == enumMapAlgorithm.end())
        {
            cout << "Algorithm: " << value << endl;
            throw runtime_error("Invalid Algorithm");
        }
        return (iValue->second);
    }

    LoadType ParseLoadType(const string& value)
    {
        EnumMapLoad::const_iterator iValue = enumMapLoad.find(value);
        if (iValue == enumMapLoad.end())
        {
            cerr << "Load Type: " << value << endl;
            throw runtime_error("Invalid Load Type");
        }
        return (iValue->second);
    }

    LoadBalancingScheme ParseLoadBalancingScheme(const string& value)
    {
        EnumMapLoadBalancingScheme::const_iterator iValue = enumMapLoadBalancingScheme.find(value);
        if (iValue == enumMapLoadBalancingScheme.end())
        {
            cerr << "Load Balancing Scheme: " << value << endl;
            throw runtime_error("Invalid Load Balancing Scheme");
        }
        return (iValue->second);
    }
};

#endif /* CONFIGURATION_HPP_ */
