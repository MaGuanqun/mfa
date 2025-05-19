#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <Eigen/Dense>

#include"ori_function.h"

namespace convert
{

    // convert matrix to vector except the last column
    template<class T>
    void convertMatrixToVect(std::vector<std::vector<VectorX<T>>>& root_vector,std::vector<MatrixX<T>>& root)
    {
        root_vector.resize(root.size());
        for(int i=0;i<root_vector.size();i++)
        {
            for(int j=0;j<root[i].rows();j++)
            {
                root_vector[i].emplace_back(root[i].block(j,0,1,root[i].cols()-1).transpose());
            }
        }

    }


    template<typename T>
    void get_derivative_value(std::vector<std::vector<T>>& derivative, std::vector<VectorX<T>>& root)
    {
        derivative.resize(root.size());
        VectorXi deriv_order(root[0].size());
        for(int i=0;i<derivative.size();++i)
        {
            derivative[i].resize(root[i].size());
            for(int j=0;j<root[i].size();++j)
            {
                deriv_order.setZero();
                deriv_order[j]=1;
                derivative[i][j]=ori_function::ackley(root[i],deriv_order);
            }
        }
    }

}