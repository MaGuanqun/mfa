#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"

#include "utility_function.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

namespace critical_point_index{

template<class T>
void compute_dev_f(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data, VectorX<T>& p,MatrixX<T>& dev_f)
{
    dev_f.resize(p.size(),p.size());
    VectorXi deriv(p.size());
    VectorX<T> dev_f_vector(p.size());    
    for(int i=0;i<p.size();i++)
    {
        for(int j=i;j<p.size();j++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[j]+=1;
            mfa->DecodePt(*mfa_data,p,deriv,dev_f_vector);
            dev_f(j,i) = dev_f_vector[0];// / (local_domain_range[i]*local_domain_range[j]);
            dev_f(i,j) = dev_f_vector[0];
        }
    }


}

template<class T>
int check_index(MatrixX<T>& matrix, T threshold)
{
    if(abs(matrix.determinant()) < threshold)
    {
        return -1;
    }

    Eigen::EigenSolver<MatrixX<T>> es(matrix);
    auto eigenvalues = es.eigenvalues();

    // std::cout<<eigenvalues[0].real()<<" "<<eigenvalues[1].real()<<std::endl;

    // Determine the type of critical point based on eigenvalues
    if(eigenvalues[0].real() > 0 && eigenvalues[1].real() > 0)
    {
        return 3;
    }
    if(eigenvalues[0].real() < 0 && eigenvalues[1].real() < 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

template<class T>
void critical_point_index(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data, std::vector<VectorX<T>>& p, std::vector<int>& index, T threshold)
{
    MatrixX<T> dev_f;
    index.resize(p.size());
    for(int i=0;i<p.size();i++)
    {        
        compute_dev_f(mfa,mfa_data,p[i],dev_f);
        index[i] = check_index(dev_f,threshold);
    }

}


void output_index_summary(std::vector<int>& index)
{
    std::vector<int> index_summary(5,0);
    for(int i=0;i<index.size();i++)
    {
        if(index[i] == -1)
        {
            index_summary[4]++;
            continue;
        }
        index_summary[index[i]]++;
    }
    std::cout << "Critical point index summary: " << std::endl;
    std::cout << "Saddle points: " << index_summary[1] << std::endl;
    std::cout << "Local minima: " << index_summary[0] << std::endl;
    std::cout << "Local maxima: " << index_summary[3] << std::endl;
    std::cout << "Degenerate: " << index_summary[4] << std::endl;

}
}