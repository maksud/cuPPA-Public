/*********************************
Developer: Maksudul Alam
Oak Ridge National Laboratory
*********************************/

#ifndef LINEARVERTEX_H_
#define LINEARVERTEX_H_

#include "Vertex.h"
#include <iostream>

#define NULL_PTR -1
#define LEFT_CHILD 0
#define RIGHT_CHILD 1

class LinearVertex: public Vertex
{
private:
    EdgeIndexType _M;
    EdgeIndexType _root;
    EdgeIndexType _size;
    EdgeIndexType _count;

    NodeType* _A;
    EdgeIndexType** _child;

public:
    LinearVertex(EdgeIndexType d);
    virtual ~LinearVertex();

    void clear()
    {
        // -1 is important
        std::fill_n(_A, _M, -1);
        std::fill_n(_child[LEFT_CHILD], _M, NULL_PTR);
        std::fill_n(_child[RIGHT_CHILD], _M, NULL_PTR);
    }

    void Adj(EdgeIndexType index, NodeType v);
    NodeType Adj(EdgeIndexType index);

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
                si = _child[LEFT_CHILD][si];
            }
            else if (x > _A[si])
            {
                si = _child[RIGHT_CHILD][si];
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
                si = _child[LEFT_CHILD][si];
            }
            else if (x > _A[si])
            {
                si = _child[RIGHT_CHILD][si];
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
                if (_child[LEFT_CHILD][ci] == NULL_PTR)
                {
                    _A[si] = v;
                    _child[LEFT_CHILD][ci] = si;

                    _count++;
                    return _count - 1;
                }
                else
                {
                    ci = _child[LEFT_CHILD][ci];
                }
            }
            else if (v > _A[ci])
            { //Go RIGHT
                if (_child[RIGHT_CHILD][ci] == NULL_PTR)
                {
                    _A[si] = v;
                    _child[RIGHT_CHILD][ci] = si;

                    _count++;
                    return _count - 1;
                }
                else
                {
                    ci = _child[RIGHT_CHILD][ci];
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
    inline bool uniqueInsert(EdgeIndexType si, NodeType x)
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
                if (_child[LEFT_CHILD][ci] == NULL_PTR)
                {

                    _A[si] = x;
                    _child[LEFT_CHILD][ci] = si;
                    return (1);
                }
                else
                {
                    ci = _child[LEFT_CHILD][ci];
                }
            }
            else if (x > _A[ci])
            { //Go RIGHT
                if (_child[RIGHT_CHILD][ci] == NULL_PTR)
                {
                    _A[si] = x;
                    _child[RIGHT_CHILD][ci] = si;
                    return (1);
                }
                else
                {
                    ci = _child[RIGHT_CHILD][ci];
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
            return removeItem(key, _child[LEFT_CHILD][node], &_child[LEFT_CHILD][node]);
        }
        else if (key > _A[node])
        { //Item in Right Sub-Tree
            return removeItem(key, _child[RIGHT_CHILD][node], &_child[RIGHT_CHILD][node]);
        }
        else
        {
            if (_child[LEFT_CHILD][node] == NULL_PTR)
            {
                *parent = _child[RIGHT_CHILD][node];
                _child[RIGHT_CHILD][node] = NULL_PTR;
                return node;
            }
            else if (_child[RIGHT_CHILD][node] == NULL_PTR)
            {
                *parent = _child[LEFT_CHILD][node];
                _child[LEFT_CHILD][node] = NULL_PTR;
                _A[node] = NULL_PTR;
                return node;
            }
            else
            { //Find Leftmost child of right child
                EdgeIndexType replacementItem = _child[RIGHT_CHILD][node];
                EdgeIndexType parentReplacementItem = NULL_PTR;
                while (_child[LEFT_CHILD][replacementItem] != NULL_PTR)
                {
                    parentReplacementItem = replacementItem;
                    replacementItem = _child[LEFT_CHILD][replacementItem];
                }

                // DOES That do all?
                *parent = replacementItem;
                _child[LEFT_CHILD][replacementItem] = _child[LEFT_CHILD][node];

                if (parentReplacementItem != NULL_PTR)
                {
                    //Find Rightmost child of replacement item
                    EdgeIndexType rightmostReplacementItem = replacementItem;
                    //EdgeIndexType parentRightmostReplacementItem = NULL_PTR;
                    while (_child[RIGHT_CHILD][rightmostReplacementItem] != NULL_PTR)
                    {
                        //parentRightmostReplacementItem = rightmostReplacementItem;
                        rightmostReplacementItem = _child[RIGHT_CHILD][rightmostReplacementItem];
                    }
                    _child[RIGHT_CHILD][rightmostReplacementItem] = _child[RIGHT_CHILD][node];
                    _child[LEFT_CHILD][parentReplacementItem] = NULL_PTR;
                }

                _child[RIGHT_CHILD][node] = NULL_PTR;
                _child[LEFT_CHILD][node] = NULL_PTR;

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
            if (uniqueInsert(index, newKey))
            {
                return index;
            }
            else
            {
                //Restore Previous
                uniqueInsert(index, oldKey);
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
            inOrderTraversal(_child[LEFT_CHILD][node]);
            std::cout << _A[node] << ", ";
            inOrderTraversal(_child[RIGHT_CHILD][node]);
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
            std::cout << i << ":\t" << _A[i] << "\t" << _child[LEFT_CHILD][i] << "\t" << _child[RIGHT_CHILD][i] << "\t" << std::endl;
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

#endif /* PBAVERTEX_H_ */
