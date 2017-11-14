/*********************************
 Developer: Maksudul Alam
 Oak Ridge National Laboratory
 *********************************/

#ifndef PAHELPER_HPP_
#define PAHELPER_HPP_

#include "PAGraph.h"

class PAHelper
{
public:

    static void outputSCMDegreeDistribution(NodeType* _G, NodeType n, EdgeIndexType m, std::string outputDir, std::string degreeFile, bool edgeList = false)
    {
        NodeType* nodeDegree = new NodeType[n];
        NodeType* degreeDist = new NodeType[n];
        for (NodeType i = 0; i < n; i++)
        {
            nodeDegree[i] = 0;
            degreeDist[i] = 0;
        }

        //Calculate Computation Degree of a Node
        for (NodeType i = m; i < n; i++)
        {
            for (EdgeIndexType j = 0; j < m; j++)
            {
                if (edgeList)
                {
                    EdgeType Z1 = i * m * 2 + j;
                    EdgeType Z2 = i * m * 2 + m + j;

                    nodeDegree[_G[Z1]] = nodeDegree[_G[Z1]] + 1;
                    nodeDegree[_G[Z2]] = nodeDegree[_G[Z2]] + 1;
                }
                else
                {
                    EdgeType Z = i * m + j;
                    nodeDegree[i] = nodeDegree[i] + 1;
                    nodeDegree[_G[Z]] = nodeDegree[_G[Z]] + 1;
                }
            }
        }

        EdgeType sum = 0;
        for (NodeType i = 0; i < n; i++)
        {
            degreeDist[nodeDegree[i]] = degreeDist[nodeDegree[i]] + 1;
            sum = sum + nodeDegree[i];
        }

        std::stringstream sstm;
        sstm << outputDir << "/" << degreeFile;
        std::ofstream outfile(sstm.str().c_str());
        outfile << "Degree,Count" << endl;
        for (NodeType i = 0; i < n; i++)
        {
            if (degreeDist[i] > 0)
                outfile << i << "," << degreeDist[i] << endl;
        }
        outfile.close();

        delete[] degreeDist;
        delete[] nodeDegree;
    }

    static void outputSBADegreeDistribution(NodeType* _G, NodeType n, EdgeIndexType m, std::string outputDir, std::string degreeFile, bool edgeList = false)
    {
        NodeType* nodeDegree = new NodeType[n];
        NodeType* degreeDist = new NodeType[n];

        for (NodeType u = 0; u < n; u++)
        {
            nodeDegree[u] = 0;
            degreeDist[u] = 0;
        }

        EdgeType maxEdges = ((EdgeType) n * m * 2);

        cout << "MAXIMUM SIZE: " << maxEdges << endl;

        //Calculate Computation Degree of a Node
        for (EdgeType u = m; u < n; u++)
        {
            for (EdgeIndexType j = 0; j < m; j++)
            {
                if (edgeList)
                {
                    EdgeType e_1 = (2 * u * m) + (2 * j);
                    EdgeType e_2 = (2 * u * m) + (2 * j) + 1;

                    if (e_1 < 0 || e_2 < 0)
                    {
                        printf("LT Zero 1\n");
                    }

                    if (e_1 >= maxEdges || e_2 >= maxEdges)
                    {
                        printf("GT Overflow 1\n");
                    }

                    if (_G[e_1] < 0 || _G[e_2] < 0)
                    {
                        printf("LT Zero 2\n");
                    }

                    if (_G[e_1] >= n || _G[e_2] >= n)
                    {
                        printf("GT Overflow 2\n");
                    }

                    nodeDegree[_G[e_1]] = nodeDegree[_G[e_1]] + 1;
                    nodeDegree[_G[e_2]] = nodeDegree[_G[e_2]] + 1;
                }
                else
                {
                    EdgeType Z = u * m + j;
                    nodeDegree[u] = nodeDegree[u] + 1;
                    nodeDegree[_G[Z]] = nodeDegree[_G[Z]] + 1;
                }
            }
        }

        printf("Degree of each vertex is computed!\n");

        EdgeType sum = 0;
        for (NodeType i = 0; i < n; i++)
        {
            degreeDist[nodeDegree[i]] = degreeDist[nodeDegree[i]] + 1;
            sum = sum + nodeDegree[i];
        }

        std::stringstream sstm;
        sstm << outputDir << "/" << degreeFile;
        std::ofstream outfile(sstm.str().c_str());
        outfile << "Degree,Count" << endl;
        for (NodeType i = 0; i < n; i++)
        {
            if (degreeDist[i] > 0)
                outfile << i << "," << degreeDist[i] << endl;
        }
        outfile.close();

        delete[] degreeDist;
        delete[] nodeDegree;
    }

    static void output_ba_graph(NodeType V[], LongLong B, ProcessorType rank, NodeType lower_bound, std::string path = "")
    {
        std::ostringstream os;
        os << path << "out.txt." << rank;

        std::ofstream outfile(os.str().c_str());
        for (LongLong i = 0; i < B; i++)
        {
            outfile << i + lower_bound << " " << V[i] << "\n";
        }
        outfile.close();
    }

    static void writeGraphToDisk(NodeType* G, NodeType N, EdgeIndexType M, std::string graphFileName)
    {
        std::ostringstream os;
        os << graphFileName;

        std::ofstream outfile(os.str().c_str());
        for (NodeType u = M; u < N; u++)
        {
            for (EdgeIndexType j = 0; j < M; j++)
            {
                outfile << u << " " << G[u * M + j] << "\n";
            }
        }
        outfile.close();
    }
    static void outputBAGraph4(PAGraph& V, LongLong N, NodeIndexType M, ProcessorType rank, NodeType lower_bound, std::string path = "")
    {
        std::ostringstream os;
        os << path << "out.txt." << rank;

        std::ofstream outfile(os.str().c_str());
        for (LongLong i = 0; i < N; i++)
        {
            if (i + lower_bound < M)
                continue;

            for (EdgeIndexType j = 0; j < M; j++)
            {
                outfile << i + lower_bound << "\t" << V.Adj(i, j) << "\n";
            }
        }
        outfile.close();
    }
    static void output_ba_graph_interleaved(NodeType** V, LongLong N, NodeIndexType M, ProcessorType P, ProcessorType rank, std::string path = "")
    {
        std::ostringstream os;
        os << path << "out.txt." << rank;

        std::ofstream outfile(os.str().c_str());
        for (LongLong i = rank; i < N; i += P)
        {
            if (i < M)
                continue;
            for (NodeIndexType j = 0; j < M; j++)
            {
                outfile << i << " " << V[i / P][j] << "\n";
            }
        }
        outfile.close();
    }

    static void print_array(std::string name, LongLong *arr, size_t size)
    {
        std::cout << name;
        for (size_t i = 0; i < size; i++)
        {
            std::cout << arr[i] << ",";
        }
        std::cout << "\n";
    }

    static void print_V(NodeType **V, LongLong N, LongLong M)
    {
        for (NodeType i = 0; i < N; i++)
        {
            for (NodeType j = 0; j < M; j++)
            {
                std::cout << i << " --> " << V[i][j] << "\t";
            }
            std::cout << std::endl;
        }
    }

    static int parseLine(char* line)
    {
        int i = strlen(line);
        while (*line < '0' || *line > '9')
            line++;
        line[i - 3] = '\0';
        i = atoi(line);
        return i;
    }

};

#endif /* PAHELPER_HPP_ */
