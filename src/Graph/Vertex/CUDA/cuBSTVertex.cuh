/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef CUBSTVERTEX_H_
#define CUBSTVERTEX_H_

#include "cuVertex.cuh"
#include <iostream>

#define NULL_PTR -1
#define LEFT_CHILD 0
#define RIGHT_CHILD 1

class cuBSTVertex: public cuVertex
{
public:
    EdgeIndexType _M;
    EdgeIndexType _count;
    EdgeIndexType _root;

    NodeType* _A;
    EdgeIndexType* _leftChild;
    EdgeIndexType* _rightChild;

public:
    CUDA_CALLABLE_MEMBER cuBSTVertex(EdgeIndexType M)
    {
        this->_M = M;
        this->_count = 0;
        this->_root = 0;
        this->_A = (NodeType*) malloc(M * sizeof(NodeType));
        this->_leftChild = (EdgeIndexType*) malloc(M * sizeof(EdgeIndexType));
        this->_rightChild = (EdgeIndexType*) malloc(M * sizeof(EdgeIndexType));

        for (int i = 0; i < _M; i++)
        {
            _A[i] = -1;
            _leftChild[i] = NULL_PTR;
            _rightChild[i] = NULL_PTR;
        }
    }

    CUDA_CALLABLE_MEMBER cuBSTVertex(EdgeIndexType M, NodeType* A)
    {
        this->_M = M;
        this->_count = 0;
        this->_root = 0;
        this->_A = A;
        this->_leftChild = (EdgeIndexType*) malloc(M * sizeof(EdgeIndexType));
        this->_rightChild = (EdgeIndexType*) malloc(M * sizeof(EdgeIndexType));
    }

    CUDA_CALLABLE_MEMBER cuBSTVertex(EdgeIndexType M, NodeType* A, EdgeIndexType* L, EdgeIndexType* R)
    {
        this->_M = M;
        this->_count = 0;
        this->_root = 0;
        this->_A = A;
        this->_leftChild = L;
        this->_rightChild = R;

        memset(_A, -1, _M * sizeof(NodeType));
        memset(_leftChild, NULL_PTR, _M * sizeof(EdgeIndexType));
        memset(_rightChild, NULL_PTR, _M * sizeof(EdgeIndexType));
    }

    CUDA_CALLABLE_MEMBER ~cuBSTVertex()
    {

    }

    CUDA_CALLABLE_MEMBER void clear()
    {
        for (int i = 0; i < _M; i++)
        {
            _A[i] = -1;
            _leftChild[i] = NULL_PTR;
            _rightChild[i] = NULL_PTR;
        }

        _root = 0;
        _count = 0;
    }

    CUDA_CALLABLE_MEMBER void Adj(EdgeIndexType index, NodeType v)
    {
        _A[index] = v;
    }

    CUDA_CALLABLE_MEMBER NodeType Adj(EdgeIndexType index)
    {
        return _A[index];
    }

    CUDA_CALLABLE_MEMBER NodeType& operator[](EdgeIndexType index)
    {
        return _A[index];
    }

    /**
     * Returns TRUE if x exists, otherwise return FALSE
     */
    CUDA_CALLABLE_MEMBER bool contains(NodeType x)
    {
        int si = _root;
        while (si >= 0)
        {
            if (x < _A[si])
            {
                si = _leftChild[si];
            }
            else if (x > _A[si])
            {
                si = _rightChild[si];
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
    CUDA_CALLABLE_MEMBER EdgeIndexType indexOf(NodeType x)
    {
        EdgeIndexType si = _root;
        while (si >= 0)
        {
            if (x < _A[si])
            {
                si = _leftChild[si];
            }
            else if (x > _A[si])
            {
                si = _rightChild[si];
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
    CUDA_CALLABLE_MEMBER bool insertItem(EdgeIndexType si, NodeType x)
    {
        if (_count == 0)
        {
            //First Element is never a duplicate, and hence the root
            _root = si;
            _A[si] = x;
            _count++;
            return true; //No Duplicate
        }

        EdgeIndexType ci = _root;
        while (ci != NULL_PTR)
        {
            if (x < _A[ci])
            {
                //Go LEFT
                if (_leftChild[ci] == NULL_PTR)
                {
                    if (_count < _M)
                    {
                        _A[si] = x;
                        _leftChild[ci] = si;
                        _count++;
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    ci = _leftChild[ci];
                }
            }
            else if (x > _A[ci])
            {
                //Go RIGHT
                if (_rightChild[ci] == NULL_PTR)
                {
                    if (_count < _M)
                    {
                        _A[si] = x;
                        _rightChild[ci] = si;
                        _count++;
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    ci = _rightChild[ci];
                }
            }
            else if (x == _A[ci])
            {
                return false;
            }
        }
        return false; //Pessimistic FAILURE
    }

    /**
     * Remove the item and returns the index of the item.
     */
    CUDA_CALLABLE_MEMBER NodeIndexType removeItem(NodeType key)
    {
        EdgeIndexType node = _root;
        EdgeIndexType *parent = &_root;

        while (node != NULL_PTR)
        {
            if (key < _A[node])
            {
                //Item in Left Sub-Tree
                parent = &_leftChild[node];
                node = _leftChild[node];
            }
            else if (key > _A[node])
            {
                //Item in Right Sub-Tree
                parent = &_rightChild[node];
                node = _rightChild[node];
            }
            else
            {
                // Found Item
                if (_leftChild[node] == NULL_PTR)
                {
                    *parent = _rightChild[node];
                    _rightChild[node] = NULL_PTR;
                    _A[node] = NULL_PTR;
                    _count--;
                    return node;
                }
                else if (_rightChild[node] == NULL_PTR)
                {
                    *parent = _leftChild[node];
                    _leftChild[node] = NULL_PTR;
                    _A[node] = NULL_PTR;
                    _count--;
                    return node;
                }
                else
                {
                    //Find Leftmost child of right child
                    EdgeIndexType replacementItem = _rightChild[node];
                    EdgeIndexType parentReplacementItem = NULL_PTR;
                    while (_leftChild[replacementItem] != NULL_PTR)
                    {
                        parentReplacementItem = replacementItem;
                        replacementItem = _leftChild[replacementItem];
                    }

                    // DOES That do all?
                    *parent = replacementItem;
                    _leftChild[replacementItem] = _leftChild[node];

                    if (parentReplacementItem != NULL_PTR)
                    {
                        //Find Rightmost child of replacement item
                        EdgeIndexType rightmostReplacementItem = replacementItem;
                        while (_rightChild[rightmostReplacementItem] != NULL_PTR)
                        {
                            rightmostReplacementItem = _rightChild[rightmostReplacementItem];
                        }
                        _rightChild[rightmostReplacementItem] = _rightChild[node];
                        _leftChild[parentReplacementItem] = NULL_PTR;
                    }

                    _A[node] = NULL_PTR;

                    _rightChild[node] = NULL_PTR;
                    _leftChild[node] = NULL_PTR;
                    _count--;

                    return node;
                }
            }
        }

        return (-1);
    }

    CUDA_CALLABLE_MEMBER EdgeIndexType updateItem(NodeType oldKey, NodeType newKey)
    {
        EdgeIndexType index = removeItem(oldKey);
        if (index != NULL_PTR)
        {
            if (insertItem(index, newKey))
            {
                return index;
            }
            else
            {
                //Restore Previous
                insertItem(index, oldKey);
                return -1;
            }
        }
        else
            return -1;
    }

    CUDA_CALLABLE_MEMBER EdgeIndexType size()
    {
        return _M;
    }

    CUDA_CALLABLE_MEMBER EdgeIndexType count()
    {
        return _count;
    }

    CUDA_CALLABLE_MEMBER void printInfo()
    {
        printf("------------------------------------------------------------------\n");
        printf("ROOT: %d\n", _root);
        printf("COUNT: %d\n", _count);
        for (EdgeIndexType i = 0; i < _M; i++)
        {
            printf("%d :\t %ld \t %d \t %d \n", i, _A[i], _leftChild[i], _rightChild[i]);
        }
        printf("------------------------------------------------------------------\n");
    }

};

#endif /* CUBSTVERTEX_H_ */
