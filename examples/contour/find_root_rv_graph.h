#pragma once

#include "mfa_extend.h"
#include "utility_function.h"
#include "ridge_valley_graph.h"


// This file is to find root of h=0 

namespace find_root_rv_graph
{    

    template <typename T>
    bool spanInRange(const mfa::MFA_Data<T>& mfa_data, VectorXi& span_index, std::vector<T>& interval)
    {
        std::vector<std::vector<T>> span_range(span_index.size());
        for(int i=0;i<span_index.size();++i)
        {    
            span_range[i].emplace_back(mfa_data.tmesh.all_knots[i][span_index[i]]);
            span_range[i].emplace_back(mfa_data.tmesh.all_knots[i][span_index[i]+1]);
        }   

        for(int i=0;i<span_index.size();++i)
        {
            if(span_range[i][0]<interval[2*i]){
                return false;
            }
            if(span_range[i][1]>interval[2*i+1]){
                return false;
            }
        }
        return true;
    }

    template<typename T>
    void compute_gradient(const Block<T>* b, VectorX<T>& p, VectorX<T>& grad)
    {
        grad.resize(p.size());
        VectorX<T> f_vector(1);
        VectorXi deriv = VectorXi::Zero(p.size());
        for(int i=0;i<p.size();i++)
        {
            deriv[i]=1;
            mfa_extend::recover_mfa(b, p,f_vector,deriv);
            grad[i] = f_vector[0];///local_domain_range[i];
            deriv[i]=0;
        }
    }

    template<typename T>
    void compute_hessian(VectorX<T>& p, const Block<T>* b,MatrixX<T>& second_deriv)
    {
        VectorXi deriv=VectorXi::Zero(p.size());
        VectorX<T> dev_f_record(b->mfa->nvars()); 
        dev_f_record.setZero();
        for(int i=0;i<p.size();i++)
        {
            deriv[i]=1;
            for(int j=i+1;j<p.size();j++)
            {
                deriv[j]=1;
                mfa_extend::recover_mfa(b, p,dev_f_record,deriv);
                second_deriv(i,j)=dev_f_record[0];
                second_deriv(j,i)=dev_f_record[0];
                deriv[j]=0;
            }
            deriv[i]=2;
            mfa_extend::recover_mfa(b,p,dev_f_record,deriv);
            second_deriv(i,i)=dev_f_record[0];
            deriv[i]=0;
        }

    }

        template<typename T>
    // only for a 2D function
    // third_deriv |f_xxx, f_xyy|
    //             |f_xxy, f_yyy| 
    void compute_third_deriv(VectorX<T>& p, const Block<T>* b, MatrixX<T>& third_deriv)
    {
        VectorXi deriv=VectorXi::Zero(p.size());
        VectorX<T> dev_f_record(b->mfa->nvars()); 
        deriv[0]=3;
        mfa_extend::recover_mfa(b,p,dev_f_record, deriv);
        third_deriv(0,0) = dev_f_record[0];

        deriv[0]=2; deriv[1]=1;
        mfa_extend::recover_mfa(b, p,dev_f_record, deriv);
        third_deriv(1,0) = dev_f_record[0];

        deriv[0]=1; deriv[1]=2;
        mfa_extend::recover_mfa(b, p, dev_f_record, deriv);
        third_deriv(0,1) = dev_f_record[0];

        deriv[0]=0; deriv[1]=3;
        mfa_extend::recover_mfa(b, p,dev_f_record, deriv);
        third_deriv(1,1) = dev_f_record[0];

    }


    template<typename T>
    // h = ||\nabla f \cross \nabla (||\nabla f||^2)
    int gradient_descent_compute_h_gradient_h(VectorX<T>& p, const Block<T>* b, VectorX<T>& h_first_deriv,T& h,T value, T threshold)
    {



        VectorX<T> f_first_deriv(p.size());
        compute_gradient(b,p,f_first_deriv);



        MatrixX<T> f_second_deriv(p.size(),p.size());
        compute_hessian(p,b,f_second_deriv);
        VectorX<T> g_first_deriv(p.size());
        ridge_valley_graph::compute_gradient_g(f_second_deriv.data(),f_first_deriv,g_first_deriv);
        h = f_first_deriv[0]*g_first_deriv[1]-f_first_deriv[1]*g_first_deriv[0];


        if(abs(h-value)<threshold)
        {
            return true;
        }

        MatrixX<T> f_third_deriv(p.size(),p.size());
        compute_third_deriv(p,b,f_third_deriv);
        MatrixX<T> g_second_deriv(p.size(),p.size());
        ridge_valley_graph::compute_hessian_g(f_second_deriv.data(),f_third_deriv.data(),f_first_deriv,g_first_deriv,g_second_deriv.data());


        h_first_deriv[0]=f_first_deriv[0]*g_second_deriv.data()[1]+g_first_deriv[1]*f_second_deriv.data()[0]-f_first_deriv[1]*g_second_deriv.data()[0]-g_first_deriv[0]*f_second_deriv.data()[1];

        h_first_deriv[1]=f_first_deriv[0]*g_second_deriv.data()[3]+g_first_deriv[1]*f_second_deriv.data()[1]-f_first_deriv[1]*g_second_deriv.data()[1]-g_first_deriv[0]*f_second_deriv.data()[3];

        return false;

    }




    template<typename T>
    //we always need to extract h=0
    bool gradient_descent_h(const Block<T>* b,VectorX<T>& result, VectorX<T>& p, int max_itr,
                std::vector<std::vector<T>>& span_range,
                T d_max_square, VectorX<T>& center,
                T root_finding_epsilon, T& function_value)
    {

        VectorX<T> grad(p.size());
        int i=0;
        while(i<max_itr)
        {

            if(gradient_descent_compute_h_gradient_h(p,b,grad,function_value,0.0,root_finding_epsilon))
            {
                if(!utility::InBlock(span_range,p))
                {
                    return false;
                }
                result = p;
                return true;
            }

            // std::cout<<"gradient descent h "<<i<<" "<<p.transpose()<<std::endl;

            T step_size=function_value/grad.squaredNorm();
            p-= step_size*grad;

            if((p-center).squaredNorm()>d_max_square)
            {
                return false;
            }
            if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
            {
                return false;
            }

            i++;

        }

        return false;

    }

    // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool root_finding(const Block<T>* b, VectorXi& span_index, std::vector<VectorX<T>>& root,
        T root_finding_epsilon, std::vector<T>& function_value, T same_root_epsilon) { 

        function_value.clear();
        root.clear();
        
        VectorXi one = VectorXi::Ones( b->mfa->var(0).p.size());
        // int deg = (mfa_data->p-one).prod();

        int maxIter=100;

        int distance_stop_itr = 5;

        std::vector<VectorX<T>> root_in_original_domain;
        
        // std::cout<<"max_iteration--"<<maxIter<<std::endl;

        // std::cout<<domain_min.transpose()<<" "<<domain_range.transpose()<<std::endl;

        std::vector<std::vector<T>> span_range(span_index.size());
        
        auto domain_range = b->core_maxs-b->core_mins;
        VectorX<T> center(span_index.size());
        for(int i=0;i<span_index.size();++i)
        {    
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]]*domain_range[i]+b->core_mins[i]);
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]+1]*domain_range[i]+b->core_mins[i]);
            center[i]=(span_range[i][0]+span_range[i][1])*0.5;
        }   


        // if(span_range[0][0]<-1.55 || span_range[0][1]>1.55 || span_range[1][0]<-0.8 || span_range[1][1]>2.3)
        // {
        //     std::cout<<"==="<<std::endl;
        //     std::cout<<span_range[0][0]<<" "<<span_range[0][1]<<" "<<span_range[1][0]<<" "<<span_range[1][1]<<std::endl;
        //     std::cout<<mfa_data->tmesh.all_knots[0][span_index[0]]<<" "<<mfa_data->tmesh.all_knots[0][span_index[0]+1]<<" "<<mfa_data->tmesh.all_knots[0][span_index[1]]<<" "<<mfa_data->tmesh.all_knots[0][span_index[1]+1]<<std::endl;
        //     std::cout<<domain_range.transpose()<<std::endl;
        //     std::cout<<domain_min.transpose()<<std::endl;
        // }
    


        // compute distance to terminate iteration
        T d_max_square=0;
        for(auto i=span_range.begin();i!=span_range.end();++i)
        {  
            d_max_square+=((*i)[1]-(*i)[0])*((*i)[1]-(*i)[0]);
        }

        d_max_square*=distance_stop_itr * distance_stop_itr; //d^2=(2*diagonal of span)^2

        std::vector<std::vector<T>>initial_point;


        utility::compute_initial_points_rv(initial_point,b->mfa->var(0).p,span_range);

        VectorXi num_initial_point_every_domain(initial_point.size());
        for(int i=0;i<num_initial_point_every_domain.size();i++)
        {
            num_initial_point_every_domain[i]=initial_point[i].size();
        }

        int num_initial_point = num_initial_point_every_domain.prod();

        VectorX<T> next_root; 




        // for(int i=0;i<initial_point.size();++i)
        // {
        //     std::cout<<initial_point[i][0]<<" "<<initial_point[i][1]<<std::endl;
        // }
        // std::cout<<std::endl;
        // std::cout << "Press Enter to continue...";
        // std::cin.get(); // Waits for the user to press Enter

        VectorXi domain_index;
        VectorXi number_in_every_domain;
        VectorX<T> current_initial_point(initial_point.size());
        utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);

        T func_value;

        for(int i=0;i<num_initial_point;++i)
        {

            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<current_initial_point.size();j++)
            {
                current_initial_point[j]=initial_point[j][domain_index[j]];
            }        



            if(gradient_descent_h(b, next_root, current_initial_point,maxIter,span_range,d_max_square,center,root_finding_epsilon,func_value))
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
                    function_value.emplace_back(func_value);
                }
                // std::cout<<"is a new root "<<std::endl;                
            }                        
        }


        return !root.empty();

    }



    template<typename T>
    bool root_finding(Block<T>* block,
    std::vector<VectorX<T>>& root,VectorXi& number_in_every_domain,
        int current_index,
        T root_finding_epsilon, std::vector<T>& function_value,T same_root_epsilon, std::vector<T>& shrink_ratio) //2^n+1 initial points) 
    {
        VectorXi span_index(block->core_maxs.size());



        for(auto i=0;i<block->mfa->nvars();++i)
        {
            utility::obtainDomainIndex(current_index,span_index,number_in_every_domain);
            span_index+=block->mfa->var(i).p;

            if(spanInRange(block->mfa->var(i),span_index,shrink_ratio))
            {
                if(root_finding(block,span_index, root,
                root_finding_epsilon,function_value,same_root_epsilon))
                {
                    return true;
                }
            }
        }

        return false;
    }


    template<typename T>
    T coeff_control_point_range(Eigen::MatrixX<T>& ctrl_pts)
    {
        // b->geometry.mfa_data->tmesh.tensor_prods[0].ctrl_pts

        // std::cout<<ctrl_pts.rows()<<" "<<ctrl_pts.cols()<<std::endl;

        T max = ctrl_pts.maxCoeff();
        T min = ctrl_pts.minCoeff();
        std::cout<<"max critical points"<<max<<std::endl;
        std::cout<<"min critical points"<<min<<std::endl;

        return 1.0/(max-min);
    }

    template<typename T>
    void span_range(Block<T>* block) 
    {
        mfa::MFA_Data<T>& mfa_data = block->mfa->var(0);
        
        VectorX<T> range(mfa_data.p.size());
        VectorX<T> range_max(mfa_data.p.size());
        range.setOnes();
        range_max.setZero();
        for(int i=mfa_data.p[0];i<mfa_data.tmesh.all_knots[0].size()-mfa_data.p[0];++i)
        {       
            for(int j=0;j<range.size();++j)
            {
                T k =  mfa_data.tmesh.all_knots[j][i+1] - mfa_data.tmesh.all_knots[j][i];
                if(k!=0.0)
                {
                    if(range[j]>k)
                        range[j]=k;
                    if(range_max[j]<k)
                        range_max[j]=k;
                } 
            }            
        }   
        std::cout<<"min span size "<<range.transpose()<<std::endl;
        std::cout<<"max span size "<<range_max.transpose()<<std::endl;
    }


}
