#pragma once

#include    <mfa/mfa.hpp>
#include    <mfa/block_base.hpp>


#include    <stdio.h>

#include    <Eigen/Dense>

#include    <random>
#include "utility_function.h"
// #include <atomic>
#include <tbb/tbb.h>

//for nd domain
namespace span_filter
{
    template <typename T>
    void valid_span(
        const Block<real_t>*              block,
        std::vector<std::vector<VectorXi>>&           valid_span,
        T func_value
    )
    {
        int nvars           = block->mfa->nvars();   
        valid_span.resize(nvars);

        for (size_t i = 0; i < nvars; i++)                                      // science variables
        {
            auto& mfa_data          = block->mfa->var(i);
            auto& tmesh             = mfa_data.tmesh;
            auto& tensor_prods      = tmesh.tensor_prods;
            auto& tc = tmesh.tensor_prods[0];
            int ndom_dims           = block->dom_dim;  // number of geometry dims
            // auto& all_knots         = tmesh.all_knots;
            // auto& all_knot_levels   = tmesh.all_knot_levels;

            // block into the original funciton,
            //the common block means the number should be the same as the domain dimension
            // memset(count_in_original_function.data(),0,sizeof(size_t)*count_in_original_function.size());

            valid_span[i].reserve(tc.nctrl_pts.prod()/10);


            int total_block=0;        

            VectorXi p =mfa_data.p;

            VectorXi control_points_num = tc.nctrl_pts;

            VectorXi span_num = control_points_num-p;

            int spanned_block_num = span_num.prod();

            std::vector<int> count_in_original_function(spanned_block_num,0);//here,we project every 

            VectorXi number_in_every_domain; //span
            utility::obtain_number_in_every_domain(span_num,number_in_every_domain);

            VectorXi block_span_num = p+Eigen::VectorXi::Ones(p.size()); //every block

            std::cout<<"block span num "<<block_span_num.transpose()<<std::endl;

            VectorXi number_in_every_dmain_block;
            utility::obtain_number_in_every_domain(block_span_num,number_in_every_dmain_block);

            VectorXi number_in_every_domain_ori; //control_point 
            utility::obtain_number_in_every_domain_row_major(control_points_num,number_in_every_domain_ori);

            //loop over every block, 
            
            int block_span=block_span_num.prod();

            // const T* control_points = mfa_data.tmesh.tensor_prods[0].ctrl_pts.data();

            tbb::affinity_partitioner ap;
            tbb::parallel_for(tbb::blocked_range<size_t>(0,spanned_block_num),
            [&](const tbb::blocked_range<size_t>& range)
            {

                for (auto j = range.begin(); j!= range.end(); j++)
                {
                    // if(j==3099){
                        VectorXi span_index(p.size());
                        //decode the span index in every dimension
                        utility::obtainDomainIndex(j,span_index,number_in_every_domain);
                        //get the index in control point for this block
                        std::vector<size_t>index_for_control_point(block_span);
                        utility::obtain_index_for_control_point(index_for_control_point,span_index,number_in_every_dmain_block,number_in_every_domain_ori);                    
                        // if(utility::check_valid_span(control_points,index_for_control_point,func_value))
                        // {
                            //project it to the original function span (fx->f)
                            // span_index+=mfa_data->p;;
                            count_in_original_function[j]=1;               
                        // }
                    // }
                }
                
            },ap
            );
            total_block = spanned_block_num;
            
            

            VectorXi index(ndom_dims);

            for(int j=0;j<count_in_original_function.size();++j)
            {
                if(count_in_original_function[j]==1)
                {
                    utility::obtainDomainIndex(j,index,number_in_every_domain);

                    valid_span[i].push_back(index);
                    

                }
                // else if(count_in_original_function[j].load()>ndom_dims){
                //     std::cout<<"error "<<count_in_original_function[j].load()<<" "<<j<<std::endl;
                // }

            }


            std::cout<<"total span num "<<total_block<<" filtered span num "<< valid_span[i].size()<<std::endl;


        }
    }



} // namespace span_filter

