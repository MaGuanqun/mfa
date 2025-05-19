#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <exception>
#include <functional>
#include <Eigen/Dense>



class KDTree
{


    public:
    template<typename T>
    KDTree(const std::vector<Eigen::VectorX<T>>& points)
    {
        std::vector<int> indices(points.size());
        for (auto i = 0; i < indices.size(); ++i) 
        {
            indices[i] = i;
        }
        root = buildTree(points,indices.data(),points.size(),0);        
    }

    template<typename T>
    void radiusSearch(const std::vector<Eigen::VectorX<T>>& points, const Eigen::VectorX<T>& query, std::vector<int>& indices, 
    T squared_R)
    {
        radiusSearchRecursive(points,query,root, indices,squared_R);
    }

    private:
    struct Node 
    {
        //Eigen::VectorX<T> point;
        int index;
        Node* left;
        Node* right;
        int axis;

        Node() : index(-1),axis(-1), left(nullptr), 
        right(nullptr) {}
    };

    Node* root;
    template<typename T>
    Node* buildTree(const std::vector<Eigen::VectorX<T>>& points, int* indices, 
    int npoints, int depth) 
    {
        if (npoints <=0 ) {
            return nullptr;
        }
        
        const int axis = depth % points[0].size(); 

        const int medianIdx = (npoints-1) / 2;
        
        std::nth_element(indices, indices + medianIdx, indices + npoints,
                            [&](int a, int b) { return points[a][axis] < points[b][axis]; });
        
        Node* node = new Node();
        node->index = indices[medianIdx];
        node->axis = axis;
        
        node->left = buildTree(points, indices, medianIdx, depth + 1);
        node->right = buildTree(points, indices+medianIdx+1, npoints - medianIdx - 1, depth + 1);
        
        return node;
    }


    void clearRecursive(Node* node)
    {
        if (node == nullptr)
            return;

        if (node->left)
            clearRecursive(node->left);

        if (node->right)
            clearRecursive(node->right);

        delete node;
    }

    template<typename T>
    void radiusSearchRecursive(const std::vector<Eigen::VectorX<T>>& points, const Eigen::VectorX<T>& query, const Node* node, std::vector<int>& indices, 
    T squared_R)
    {
        if (node == nullptr)
            return;

       // const PointT& train = points_[node->index];
        T dist =  (points[node->index] - query).squaredNorm();

        if (dist <= squared_R)
            indices.push_back(node->index);

        //const int axis = node->axis;
        int dir = query[node->axis] < points[node->index][node->axis]?0:1;

        if(dir)
        {
             radiusSearchRecursive(points, query, node->right, indices, squared_R);
        }
        else
        {
            radiusSearchRecursive(points, query, node->left, indices, squared_R);
        }        

        T diff = query[node->axis] - points[node->index][node->axis];
        if (diff*diff <= squared_R)
        {
            if(dir)
            {
                radiusSearchRecursive(points, query, node->left, indices, squared_R);
            }
            else
            {
                radiusSearchRecursive(points, query, node->right, indices, squared_R);
            }      

        }
    }
};