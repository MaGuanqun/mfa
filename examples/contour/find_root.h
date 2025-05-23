#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"
#include "utility_function.h"
#include "mfa_extend.h"

namespace find_root
{    

    template<typename T>
    void compute_gradient(const Block<T>* b, VectorX<T>& p, VectorX<T>& grad)
    {
        grad.resize(p.size());
        VectorX<T> f_vector(1);
        VectorXi deriv = VectorXi::Zero(p.size());
        for(int i=0;i<p.size();i++)
        {
            deriv[i]=1;
            mfa_extend::recover_mfa(b,p,f_vector,deriv);
            // mfa_extend::recover_mfa_selected(mfa,mfa_data,span_index, p,deriv,weights,f_vector,domain_min,domain_range);
            grad[i] = f_vector[0];///local_domain_range[i];
            deriv[i]=0;
        }
    }

    template<typename T>
    //point_func_value: the value of the point we want to extract
    bool gradient_descent(const Block<T>* b,VectorX<T>& result, VectorX<T>& p, int max_itr,
                std::vector<std::vector<T>>& span_range,
                T d_max_square, VectorX<T>& center,
                T root_finding_epsilon, VectorX<T>& function_value, T point_func_value)
    {

        mfa_extend::recover_mfa(b,p,function_value);
        function_value[0]-=point_func_value;
        if(std::abs(function_value[0])<root_finding_epsilon)
        {
            result = p;
            return true;
        }

        VectorX<T> grad;
        int i=0;
        while(i<max_itr)
        {
            compute_gradient(b, p, grad);
            T step_size=function_value[0]/grad.squaredNorm();
            // printf("step_size %f\n",step_size);
            p-= step_size*grad;
            if((p-center).squaredNorm()>d_max_square)
            {
                return false;
            }
            if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
            {
                return false;
            }
            mfa_extend::recover_mfa(b,p,function_value);
            function_value[0]-=point_func_value;
            

            if(std::abs(function_value[0])<root_finding_epsilon)
            {
                if(!utility::InBlock(span_range,p))
                {
                    return false;
                }
                result = p;
                function_value[0]+=point_func_value;
                return true;
            }
            i++;

        }

        return false;

    }

    // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool root_finding(const Block<T>* b, VectorXi& span_index, std::vector<VectorX<T>>& root,
        T root_finding_grad_epsilon, std::vector<T>& function_value, T point_func_value, T same_root_epsilon) { 

        function_value.clear();
        root.clear();
        
        VectorXi one = VectorXi::Ones(b->mfa->var(0).p.size());
        // int deg = (mfa_data->p-one).prod();

        int maxIter=30;

        int distance_stop_itr = 5;

        std::vector<VectorX<T>> root_in_original_domain;
        
        // std::cout<<"max_iteration--"<<maxIter<<std::endl;

        std::vector<std::vector<T>> span_range(span_index.size());
        
        auto domain_range = b->core_maxs-b->core_mins;

        VectorX<T> center(span_index.size());
        for(int i=0;i<span_index.size();++i)
        {    
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]]*domain_range[i]+b->core_mins[i]);
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]+1]*domain_range[i]+b->core_mins[i]);

            center[i]=(span_range[i][0]+span_range[i][1])*0.5;
        }   



        // compute distance to terminate iteration
        T d_max_square=0;
        for(auto i=span_range.begin();i!=span_range.end();++i)
        {  
            d_max_square+=((*i)[1]-(*i)[0])*((*i)[1]-(*i)[0]);
        }

        d_max_square*=distance_stop_itr * distance_stop_itr; //d^2=(2*diagonal of span)^2

        std::vector<std::vector<T>>initial_point;


        utility::compute_initial_points_iso(initial_point,b->mfa->var(0).p,span_range);

        VectorXi num_initial_point_every_domain(initial_point.size());
        for(int i=0;i<num_initial_point_every_domain.size();i++)
        {
            num_initial_point_every_domain[i]=initial_point[i].size();
        }

        int num_initial_point = num_initial_point_every_domain.prod();

        VectorX<T> next_root; 

        // std::cout<<"initial_point "<<num_initial_point<<std::endl;

        // for(int i=0;i<initial_point.size();++i)
        // {
        //     std::cout<<initial_point[i]<<" ";
        // }
        // std::cout<<std::endl;
        // std::cout << "Press Enter to continue...";
        // std::cin.get(); // Waits for the user to press Enter

        VectorXi domain_index;
        VectorXi number_in_every_domain;
        VectorX<T> current_initial_point(initial_point.size());
        utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);

        VectorX<T> func_value;

        for(int i=0;i<num_initial_point;++i)
        {

            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<current_initial_point.size();j++)
            {
                current_initial_point[j]=initial_point[j][domain_index[j]];
            }        
            // current_initial_point=ini_p;


            // std::cout<<"intial point "<< current_initial_point.transpose()<<std::endl;
            // std::cout<< "initial_point "<<i<<" "<<  current_initial_point.transpose()<<std::endl;

            if(gradient_descent(b, next_root, current_initial_point,maxIter,span_range,d_max_square,center,root_finding_grad_epsilon,func_value,point_func_value))
            {
                bool duplicate = false;
                for(auto i=root.begin();i!=root.end();++i)
                {
                    if(((*i)-next_root).squaredNorm()<same_root_epsilon*same_root_epsilon)
                    {
                        duplicate = true;
                        break;
                    }
                }
                if(!duplicate)
                {
                    root.emplace_back(next_root);
                    function_value.emplace_back(func_value[0]);
                }
                // std::cout<<"is a new root "<<std::endl;                
            }                        
        }


        return !root.empty();

    }



    template<typename T>
    bool root_finding(const Block<real_t>* block, std::vector<std::vector<VectorXi>>& span_index, 
    std::vector<VectorX<T>>& root,//std::vector<int>& multi_of_root,
        //MatrixX<T>&             ctrl_pts,   //control points of first derivative
        int current_index,
        T root_finding_epsilon, std::vector<T>& function_value, T point_func_value,T same_root_epsilon) //2^n+1 initial points) 
    {

        for(auto i=0;i<block->mfa->nvars();++i)
        {
            VectorXi span_index_local = span_index[i][current_index]+block->mfa->var(0).p;

            if(root_finding(block,span_index_local, root,
            root_finding_epsilon,function_value,point_func_value,same_root_epsilon))
            {
                return true;
            }
        }

        return false;
    }
}
