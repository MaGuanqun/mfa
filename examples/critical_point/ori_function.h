#pragma once

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include    <stdio.h>
#include    <random>
#include "block.hpp"
#include <Eigen/Dense>


namespace ori_function
{
//first_derivative_index: -1 origin function, 0 partial x, 1 partial y
template <typename T>
T sinc(VectorX<T>& x, VectorXi& derivative_order)
{   
    T retval = 1.0;
    int dom_dim=x.size(); 
    for(int i=0;i<dom_dim;i++)
    {
        if (x[i] != 0.0)
        {

            switch (derivative_order[i])
            {
            case 0:
                retval *=sin(x[i]) / x[i];
                break;
            
            case 1:
                retval*= (cos(x[i])/x[i]-sin(x[i])/(x[i]*x[i]));//the first derivative
                break;
            
            case 2:
                retval*=(sin(x[i])*(2.0/(x[i]*x[i]*x[i])-1.0/x[i])-2.0*cos(x[i])/(x[i]*x[i]));        
                break;
            }

        }
        else
        {
            switch (derivative_order[i])
            {
            case 1:
                retval*=0.0;
                break;
            case 2:
                retval *=-1.0/3.0;
                break;
            }
        }

    }
    return 10.0*retval;

}


template <typename T>
T polysinc1(VectorX<T>& domain_pt, VectorXi& derivative_order)
{
     // a = (x + 1)^2 + (y - 1)^2 + (z + 1)^2 + ...
        // b = (x - 1)^2 + (y + 1)^2 + (z - 1)^2 + ...
    // a1 = sinc(a); b1 = sinc(b)

    T a = 0.0;
    T b = 0.0;
    for (auto i = 0; i < domain_pt.size(); i++)
    {
        T s, r;
        if (i % 2 == 0)
        {
            s = domain_pt(i) + 1.0;
            r = domain_pt(i) - 1.0;
        }
        else
        {
            s = domain_pt(i) - 1.0;
            r = domain_pt(i) + 1.0;
        }
        a += (s * s);
        b += (r * r);
    }

    return 1.0;
}

template <typename T>
T ackley(VectorX<T>& domain_pt, VectorXi& derivative_order)
{
    T a=20;
    T b=0.2;
    T c=2.0*M_PI;
    T d=domain_pt.size();

    
    T term_0=-b*sqrt(domain_pt.squaredNorm()/d);
    T term_1=0.0;
    for(int i=0;i<domain_pt.size();i++)
    {
        term_1+=cos(c*domain_pt[i]);
    }
    term_1/=d;
    
    term_0 = exp(term_0);
    term_1 = exp(term_1);

    T sum= derivative_order.sum();
    if(sum==0)
    {
       return -a*term_0-term_1+a+M_E;
    }
    else if(sum==1)
    {
        int derivative_index=0;
        for(int i=0;i<derivative_order.size();i++)
        {
            if(derivative_order[i]==1)
            {
                derivative_index = i;
                break;
            }
        }

        T partial_0=-b*domain_pt[derivative_index]/sqrt(d*domain_pt.squaredNorm());
        T partial_1=-c/d*sin(c*domain_pt[derivative_index]);

        return  -a*term_0*partial_0-term_1*partial_1;

    }
    else
    {
        int derivative_index=-1;
        for(int i=0;i<derivative_order.size();i++)
        {
            if(derivative_order[i]==2)
            {
                derivative_index = i;
                break;
            }
        }

        T prod = domain_pt.norm();

        if(derivative_index!=-1)
        {
            T partial_0=-b*domain_pt[derivative_index]/sqrt(d*domain_pt.squaredNorm());
            T partial_1=-c/d*sin(c*domain_pt[derivative_index]);
            T partial_0_2 = -b/sqrt(d)*(1.0/prod-domain_pt[derivative_index]*domain_pt[derivative_index]/(prod*prod*prod));
            T partial_1_2 = -c*c/d*cos(c*domain_pt[derivative_index]);
            return -a*term_0*(partial_0*partial_0+partial_0_2)-term_1*(partial_1*partial_1+partial_1_2);           
        }
        
        std::vector<int> derive_index_2;
        for(int i=0;i<derivative_order.size();i++)
        {
            if(derivative_order[i]==1)
            {
                derive_index_2.emplace_back(i);
            }
        }

        T partial_0_0=-b*domain_pt[derive_index_2[0]]/sqrt(d*domain_pt.squaredNorm());
        T partial_0_1=-b*domain_pt[derive_index_2[1]]/sqrt(d*domain_pt.squaredNorm());

        T partial_1_0=-c/d*sin(c*domain_pt[derive_index_2[0]]);
        T partial_1_1=-c/d*sin(c*domain_pt[derive_index_2[1]]);

        T partial_0_2 = b/sqrt(d)*partial_0_0*partial_0_1/(prod*prod*prod);
        T partial_1_2 = 0.0;

        return -a*term_0*(partial_0_0*partial_0_1+partial_0_2)-term_1*partial_1_0*partial_1_1;

    }
    
    return 0;

}


}

