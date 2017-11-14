/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef BSTVERTEX_H_
#define BSTVERTEX_H_

#include "Vertex.h"
#include <iostream>

#define NULL_PTR -1
#define LEFT_CHILD 0
#define RIGHT_CHILD 1

class BSTVertex: public Vertex
{
private:
    EdgeIndexType _M;
    EdgeIndexType _root;
    EdgeIndexType _size;
    EdgeIndexType _count;

    NodeType* _A;

    EdgeIndexType* _leftChild;
    EdgeIndexType* _rightChild;

public:

    BSTVertex(EdgeIndexType M)
    {
        this->_root = 0;
        this->_M = M;
        this->_size = M;
        this->_count = 0;
        this->_A = new NodeType[M];
        this->_leftChild = (EdgeIndexType*) malloc(M * sizeof(EdgeIndexType));
        this->_rightChild = (EdgeIndexType*) malloc(M * sizeof(EdgeIndexType));

        // -1 is important
        std::fill_n(_A, _M, -1);
        std::fill_n(_leftChild, _M, NULL_PTR);
        std::fill_n(_rightChild, _M, NULL_PTR);
    }

    virtual ~BSTVertex()
    {
        delete[] _leftChild;
        delete[] _rightChild;
        delete[] _A;
    }

    void clear()
    {
        // -1 is important
        std::fill_n(_A, _M, -1);
        std::fill_n(_leftChild, _M, NULL_PTR);
        std::fill_n(_rightChild, _M, NULL_PTR);
        _root = 0;
        _count = 0;
    }

    void Adj(EdgeIndexType index, NodeType v)
    {
        _A[index] = v;
    }

    NodeType Adj(EdgeIndexType index)
    {
        return _A[index];
    }

    NodeType& operator[](EdgeIndexType index)
    {
        return _A[index];
    }

    /**
     * Returns TRUE if x exists, otherwise return FALSE
     */
    bool contains(NodeType x)
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
    EdgeIndexType indexOf(NodeType x)
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

    inline void insert(EdgeIndexType si, NodeType x)
    {
        _A[si] = x;
    }

    EdgeIndexType insertItem(NodeType v)
    {
        EdgeIndexType si = _count;

        if (si == _root)
        { //First Element is never a duplicate, and hence the root
            _root = si;
            _A[si] = v;

            _count++;
            return _count - 1;
        }

        EdgeIndexType ci = _root;
        while (ci != NULL_PTR)
        {
            if (v < _A[ci])
            { //Go LEFT
                if (_leftChild[ci] == NULL_PTR)
                {
                    _A[si] = v;
                    _leftChild[ci] = si;

                    _count++;
                    return _count - 1;
                }
                else
                {
                    ci = _leftChild[ci];
                }
            }
            else if (v > _A[ci])
            { //Go RIGHT
                if (_rightChild[ci] == NULL_PTR)
                {
                    _A[si] = v;
                    _rightChild[ci] = si;

                    _count++;
                    return _count - 1;
                }
                else
                {
                    ci = _rightChild[ci];
                }
            }
            else if (v == _A[ci])
            {
                return -1;
            }
        }
        return -1; //Pessimistic FAILURE
    }

    /**
     * Returns TRUE if uniqueInsert is SUCCESS, otherwise return FALSE
     */
    inline bool insertItem(EdgeIndexType si, NodeType x)
    {
        if (si == _root)
        { //First Element is never a duplicate, and hence the root
            _root = si;
            _A[si] = x;
            return (1); //No Duplicate
        }

        EdgeIndexType ci = _root;
        while (ci != NULL_PTR)
        {
            if (x < _A[ci])
            { //Go LEFT
                if (_leftChild[ci] == NULL_PTR)
                {

                    _A[si] = x;
                    _leftChild[ci] = si;
                    return (1);
                }
                else
                {
                    ci = _leftChild[ci];
                }
            }
            else if (x > _A[ci])
            { //Go RIGHT
                if (_rightChild[ci] == NULL_PTR)
                {
                    _A[si] = x;
                    _rightChild[ci] = si;
                    return (1);
                }
                else
                {
                    ci = _rightChild[ci];
                }
            }
            else if (x == _A[ci])
            { //Duplicate
//                cout << "Duplicate";
                return (0);
            }
        }
        return (0); //Pessimistic FAILURE
    }

    /**
     * Returns index of the Item removed
     */
    EdgeIndexType removeItem(NodeType key, EdgeIndexType node, EdgeIndexType* parent)
    {
        if (node == NULL_PTR)
        { //Item not in BST
            return (-1);
        }

        if (key < _A[node])
        { //Item in Left Sub-Tree
            return removeItem(key, _leftChild[node], &_leftChild[node]);
        }
        else if (key > _A[node])
        { //Item in Right Sub-Tree
            return removeItem(key, _rightChild[node], &_rightChild[node]);
        }
        else
        {
            if (_leftChild[node] == NULL_PTR)
            {
                *parent = _rightChild[node];
                _rightChild[node] = NULL_PTR;
                return node;
            }
            else if (_rightChild[node] == NULL_PTR)
            {
                *parent = _leftChild[node];
                _leftChild[node] = NULL_PTR;
                _A[node] = NULL_PTR;
                return node;
            }
            else
            { //Find Leftmost child of right child
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
                    //EdgeIndexType parentRightmostReplacementItem = NULL_PTR;
                    while (_rightChild[rightmostReplacementItem] != NULL_PTR)
                    {
                        //parentRightmostReplacementItem = rightmostReplacementItem;
                        rightmostReplacementItem = _rightChild[rightmostReplacementItem];
                    }
                    _rightChild[rightmostReplacementItem] = _rightChild[node];
                    _leftChild[parentReplacementItem] = NULL_PTR;
                }

                _rightChild[node] = NULL_PTR;
                _leftChild[node] = NULL_PTR;

                return node;
            }
        }
        return (-1);
    }

    EdgeIndexType removeItem(NodeType key)
    {
        return removeItem(key, _root, &_root);
    }

    EdgeIndexType updateItem(NodeType oldKey, NodeType newKey)
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

    void inOrderTraversal(EdgeIndexType node)
    {
        if (node >= 0)
        {
            inOrderTraversal(_leftChild[node]);
            std::cout << _A[node] << ", ";
            inOrderTraversal(_rightChild[node]);
        }
    }
    void inOrderTraversal()
    {
        inOrderTraversal(_root);
    }

    void printInfo()
    {
        std::cout << std::endl << "------------------------------------------------------------------" << std::endl;
        std::cout << "ROOT: " << _root << std::endl;
        for (EdgeIndexType i = 0; i < _M; i++)
        {
            std::cout << i << ":\t" << _A[i] << "\t" << _leftChild[i] << "\t" << _rightChild[i] << "\t" << std::endl;
        }
        std::cout << std::endl << "------------------------------------------------------------------" << std::endl;
        inOrderTraversal();
    }

    EdgeIndexType size()
    {
        return _size;
    }

    EdgeIndexType count()
    {
        return _count;
    }

};

#endif /* BSTVERTEX_H_ */
