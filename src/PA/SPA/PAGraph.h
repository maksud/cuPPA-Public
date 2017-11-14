/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef GRAPH_H_
#define GRAPH_H_

#include "../../Common/DataTypes.h"
#include <iostream>
#include <climits>
using namespace std;

#define NULL_PTR -1
#define LEFT_CHILD 0
#define RIGHT_CHILD 1

class PAGraph
{
private:
    NodeIndexType _N;
    EdgeIndexType _M;

    EdgeIndexType* _root;
    EdgeIndexType* _size;
    NodeType** _V;
    EdgeIndexType*** _child;

public:
    PAGraph(NodeIndexType n, EdgeIndexType d);
    virtual ~PAGraph();

    void Adj(NodeIndexType node, EdgeIndexType edge, NodeType v);
    NodeType Adj(NodeIndexType node, EdgeIndexType edge);

    NodeType*& operator[](NodeIndexType index)
    {
        return _V[index];
    }

    /**
     * Returns TRUE if x exists, otherwise return FALSE
     */
    bool contains(NodeIndexType nodeIndex, NodeType x)
    {
        EdgeIndexType si = _root[nodeIndex];
        while (si >= 0)
        {
            if (x < _V[nodeIndex][si])
            {
                si = _child[nodeIndex][LEFT_CHILD][si];
            }
            else if (x > _V[nodeIndex][si])
            {
                si = _child[nodeIndex][RIGHT_CHILD][si];
            }
            else
            {
                return (1); //Found Item
            }
        }
        return (0);
    }

    /**
     * Return index if x exists, otherwise return -1
     */
    EdgeIndexType searchItem(NodeIndexType nodeIndex, NodeType x)
    {
        EdgeIndexType si = _root[nodeIndex];
        while (si >= 0)
        {
            if (x < _V[nodeIndex][si])
            {
                si = _child[nodeIndex][LEFT_CHILD][si];
            }
            else if (x > _V[nodeIndex][si])
            {
                si = _child[nodeIndex][RIGHT_CHILD][si];
            }
            else
            {
                return si; //Found Item
            }
        }
        return -1;
    }

    /**
     * Returns TRUE if uniqueInsert is SUCCESS, otherwise return FALSE
     */
    bool uniqueInsert(NodeIndexType nodeIndex, EdgeIndexType si, NodeType x)
    {
        if (_size[nodeIndex] == 0)
        {
            //First Element is never a duplicate, and hence the root
            _root[nodeIndex] = si;
            _V[nodeIndex][si] = x;
            _size[nodeIndex]++;
            return (1); //No Duplicate
        }

        EdgeIndexType ci = _root[nodeIndex];
        while (ci != NULL_PTR)
        {
            if (x < _V[nodeIndex][ci]) //Go LEFT
            {
                if (_child[nodeIndex][LEFT_CHILD][ci] == NULL_PTR)
                {
                    if (_size[nodeIndex] < _M)
                    {
                        _V[nodeIndex][si] = x;
                        _child[nodeIndex][LEFT_CHILD][ci] = si;
                        _size[nodeIndex]++;
                        return (1);
                    }
                    else
                    {
                        return (0);
                    }
                }
                else
                    ci = _child[nodeIndex][LEFT_CHILD][ci];
            }
            else if (x > _V[nodeIndex][ci]) //Go RIGHT
            {
                if (_child[nodeIndex][RIGHT_CHILD][ci] == NULL_PTR)
                {
                    if (_size[nodeIndex] < _M)
                    {
                        _V[nodeIndex][si] = x;
                        _child[nodeIndex][RIGHT_CHILD][ci] = si;
                        _size[nodeIndex]++;
                        return (1);
                    }
                    else
                        return (0);
                }
                else
                    ci = _child[nodeIndex][RIGHT_CHILD][ci];
            }
            else if (x == _V[nodeIndex][ci]) //Duplicate
            {
                return (0);
            }
        }
        return (0); //Pessimistic FAILURE
    }

    /**
     * Returns index of the Item removed
     */
    EdgeIndexType removeItem(NodeIndexType nodeIndex, NodeType key, EdgeIndexType ptr, EdgeIndexType* parent)
    {
        if (ptr == NULL_PTR) //Item not in BST
        {
            return (-1);
        }

        if (key < _V[nodeIndex][ptr]) //Item in Left Sub-Tree
        {
            return removeItem(nodeIndex, key, _child[nodeIndex][LEFT_CHILD][ptr], &_child[nodeIndex][LEFT_CHILD][ptr]);
        }
        else if (key > _V[nodeIndex][ptr]) //Item in Right Sub-Tree
        {
            return removeItem(nodeIndex, key, _child[nodeIndex][RIGHT_CHILD][ptr], &_child[nodeIndex][RIGHT_CHILD][ptr]);
        }
        else
        {
            if (_child[nodeIndex][LEFT_CHILD][ptr] == NULL_PTR)
            {
                *parent = _child[nodeIndex][RIGHT_CHILD][ptr];
                _child[nodeIndex][RIGHT_CHILD][ptr] = NULL_PTR;
                _size[nodeIndex]--;
                cout << "EMPTY INDEX: " << ptr << endl;
                return ptr;
            }
            else if (_child[nodeIndex][RIGHT_CHILD][ptr] == NULL_PTR)
            {
                *parent = _child[nodeIndex][LEFT_CHILD][ptr];
                _child[nodeIndex][LEFT_CHILD][ptr] = NULL_PTR;
                _size[nodeIndex]--;
                cout << "EMPTY INDEX: " << ptr << endl;
                return ptr;
            }
            else
            {
                cout << "@@@@@@@@@@@@@@@@@" << endl;
                //Find Left Most Child of Right-Child
                EdgeIndexType tmp = _child[nodeIndex][RIGHT_CHILD][ptr];
                EdgeIndexType par = NULL_PTR;
                while (_child[nodeIndex][LEFT_CHILD][tmp] != NULL_PTR)
                {
                    par = tmp;
                    tmp = _child[nodeIndex][LEFT_CHILD][tmp];
                }
                _V[nodeIndex][ptr] = _V[nodeIndex][tmp];
                if (par != NULL_PTR)
                {
                    //tmp is a LEFT child of it's parent
                    return removeItem(nodeIndex, _V[nodeIndex][_child[nodeIndex][LEFT_CHILD][par]], _child[nodeIndex][LEFT_CHILD][par], &_child[nodeIndex][LEFT_CHILD][par]);
                }
                else
                {
                    //tmp is NOT a LEFT child of it's parent
                    return removeItem(nodeIndex, _V[nodeIndex][_child[nodeIndex][RIGHT_CHILD][ptr]], _child[nodeIndex][RIGHT_CHILD][ptr], &_child[nodeIndex][RIGHT_CHILD][ptr]);
                }
            }

        }
        return (-1);
    }

    EdgeIndexType removeItem(NodeIndexType nodeIndex, NodeType key)
    {
        return removeItem(nodeIndex, key, _root[nodeIndex], &_root[nodeIndex]);
    }

    EdgeIndexType updateItem(NodeIndexType nodeIndex, NodeType oldKey, NodeType newKey)
    {
        EdgeIndexType index = removeItem(nodeIndex, oldKey);
        if (index != NULL_PTR)
        {
            if (uniqueInsert(nodeIndex, index, newKey))
            {
                return index;
            }
            else
            {
                //Restore Previous
                uniqueInsert(nodeIndex, index, oldKey);
                return -1;
            }
        }
        else
            return -1;
    }

    void inOrderTraversal(NodeIndexType nodeIndex, EdgeIndexType node)
    {
        if (node >= 0)
        {
            inOrderTraversal(nodeIndex, _child[nodeIndex][LEFT_CHILD][node]);
            cout << _V[nodeIndex][node] << ", ";
            inOrderTraversal(nodeIndex, _child[nodeIndex][RIGHT_CHILD][node]);
        }
    }

    void inOrderTraversal(NodeIndexType nodeIndex)
    {
        if (_size[nodeIndex] > 0)
        {
            inOrderTraversal(nodeIndex, _root[nodeIndex]);
        }
    }

    void printInfo(NodeIndexType nodeIndex)
    {
        cout << endl << "------------------------------------------------------------------" << endl;
        cout << "ROOT: " << _root[nodeIndex] << "\tSIZE: " << _size[nodeIndex] << endl;
        for (EdgeIndexType i = 0; i < _M; i++)
        {
            cout << i << ":\t" << _V[nodeIndex][i] << "\t" << _child[nodeIndex][LEFT_CHILD][i] << "\t" << _child[nodeIndex][RIGHT_CHILD][i] << "\t" << endl;
        }
        cout << endl << "------------------------------------------------------------------" << endl;
        inOrderTraversal(nodeIndex);
        cout << endl;
    }
};

#endif /* GRAPH_H_ */
