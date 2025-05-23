#pragma once

// #include <iostream>
// #include <complex>
// #include <vector>
// #include <cmath>
// #include <map>

// #include <mfa/mfa.hpp>



// #include "opts.h"

// #include "block.hpp"


#include "mfa_extend.h"

#include "utility_function.h"
#include "find_same_center.h"
#include "transfer_data.h"
#include "trace_in_span.h"

#include "trace_in_span.h"

// template<typename T>
// struct TraceInSpan
// {
//     size_t span_index;
//     std::vector<std::vector<std::vector<VectorX<T>>>> trace;
//     // [different values][different traces][points in a trace]
//     std::vector<std::vector<std::vector<T>>> gradient_magnitude;
//     // corresponding gradient magnitude
//     std::vector<std::vector<bool>>is_loop;
//     // whether the trace in this span is already a loop or not [different values][different traces]
//     std::vector<T> trace_value;
//     // the corresponding value
//     std::vector<std::vector<std::vector<T>>> actual_trace_value;
//     // the actual value of every point in the trace
// };

namespace jacobi_set
{

    template<typename T>
    void getFunctionValue(const Block<T>* b, std::vector<std::vector<VectorX<T>>>& root, std::vector<std::vector<T>>& value)
    {
        value.resize(root.size());
        VectorX<T> value0(1); 
        for(auto i=0;i<root.size();++i)
        {
            value[i].resize(root[i].size());
            for(auto j=0;j<root[i].size();++j)
            {
                mfa_extend::recover_mfa(b, root[i][j],value0);
                // VectorX<T> param = (root[i][j]-domain_min).cwiseQuotient(domain_range);
                // mfa->DecodePt(*mfa_data,param,value0);
                value[i][j]=value0[0];
            }
        }      
    
    }

    template<typename T>
    void getFunctionValue(const Block<T>* b, std::vector<TraceInSpan<T>>& traces_in_span)
    {
        tbb::affinity_partitioner ap;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, traces_in_span.size()),
            [&](const tbb::blocked_range<size_t>& range)
            {
                for(auto i=range.begin();i!=range.end();++i)
                {
                    getFunctionValue(b,traces_in_span[i].traces,traces_in_span[i].traces_values);
                }
            },ap
        );

    }




    template<typename T>
    void compute_gradient(VectorX<T>& p, const Block<T>* b, VectorX<T>& first_deriv)
    {
        VectorXi deriv=VectorXi::Zero(p.size());
        VectorX<T> dev_f_record(1); 
        for(int i=0;i<p.size();i++)
        {
            deriv[i]=1;
            mfa_extend::recover_mfa(b, p,dev_f_record,deriv);
            // std::cout<<"dev_f_record(0) "<<deriv.transpose()<<" | "<<p.traspose()<<" | "<< dev_f_record.transpose()<<std::endl;
            // pause();
            first_deriv[i]=dev_f_record[0];
            deriv[i]=0;
        }
    }


    template<typename T>
    void compute_hessian(VectorX<T>& p, const Block<T>* b, MatrixX<T>& second_deriv)
    {
        VectorXi deriv=VectorXi::Zero(p.size());
        VectorX<T> dev_f_record(1); 
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
            mfa_extend::recover_mfa(b, p,dev_f_record,deriv);
            second_deriv(i,i)=dev_f_record[0];
            deriv[i]=0;
        }

    }

    template<typename T>
    void compute_deriv_magnitude(VectorX<T>& p, const Block<T>* b,T& gradient_mag, T& second_deriv_tangent)
    {
        VectorX<T> first_deriv(p.size());
        compute_gradient(p,b,first_deriv);
        MatrixX<T> second_deriv(p.size(),p.size());
        compute_hessian(p,b,second_deriv);
        gradient_mag = first_deriv.squaredNorm();
        VectorX<T> tangent_direction(p.size());
        tangent_direction[0]=first_deriv[1];
        tangent_direction[1]=-first_deriv[0];
        // second_deriv_tangent = tangent_direction.dot(second_deriv*tangent_direction);

    }


    template<typename T>
    void compute_gradient_h(VectorX<T>& p, const Block<T>* b, 
    const Block<T>* b2,
    VectorX<T>& h_first_deriv, T& gradient_magnitude, T& second_deriv_tangent,
    bool need_magnitude=false)
    {
        VectorX<T> f_first_deriv(p.size());
        compute_gradient(p,b,f_first_deriv);
        MatrixX<T> f_second_deriv(p.size(),p.size());
        compute_hessian(p,b,f_second_deriv);
        

        VectorX<T> g_first_deriv(p.size());
        compute_gradient(p,b2,g_first_deriv);
        MatrixX<T> g_second_deriv(p.size(),p.size());
        compute_hessian(p,b2,g_second_deriv);


        h_first_deriv[0]=f_first_deriv[0]*g_second_deriv.data()[1]+g_first_deriv[1]*f_second_deriv.data()[0]-f_first_deriv[1]*g_second_deriv.data()[0]-g_first_deriv[0]*f_second_deriv.data()[1];

        h_first_deriv[1]=f_first_deriv[0]*g_second_deriv.data()[3]+g_first_deriv[1]*f_second_deriv.data()[1]-f_first_deriv[1]*g_second_deriv.data()[1]-g_first_deriv[0]*f_second_deriv.data()[3];

    }



    template<typename T>
    // h = ||\nabla f \cross \nabla (||\nabla f||^2)
    void compute_h(VectorX<T>& p,const Block<T>* b, const Block<T>* b2, T& h)
    {
        VectorX<T> f_first_deriv(p.size());
        compute_gradient(p,b,f_first_deriv);
        VectorX<T> g_first_deriv(p.size());
        compute_gradient(p,b2,g_first_deriv);

        h = f_first_deriv[0]*g_first_deriv[1]-f_first_deriv[1]*g_first_deriv[0];

    }

    template<typename T>
    void compute_h_gradient_h(VectorX<T>& p, const Block<T>* b,
    const Block<T>* b2, T& h, VectorX<T>& h_first_deriv)
    {
        VectorX<T> f_first_deriv(p.size());
        compute_gradient(p,b,f_first_deriv);
        MatrixX<T> f_second_deriv(p.size(),p.size());
        compute_hessian(p,b,f_second_deriv);

        VectorX<T> g_first_deriv(p.size());
        compute_gradient(p,b2,g_first_deriv);
        MatrixX<T> g_second_deriv(p.size(),p.size());
        compute_hessian(p,b2,g_second_deriv);

        h = f_first_deriv[0]*g_first_deriv[1]-f_first_deriv[1]*g_first_deriv[0];

        h_first_deriv[0]=f_first_deriv[0]*g_second_deriv.data()[1]+g_first_deriv[1]*f_second_deriv.data()[0]-f_first_deriv[1]*g_second_deriv.data()[0]-g_first_deriv[0]*f_second_deriv.data()[1];

        h_first_deriv[1]=f_first_deriv[0]*g_second_deriv.data()[3]+g_first_deriv[1]*f_second_deriv.data()[1]-f_first_deriv[1]*g_second_deriv.data()[1]-g_first_deriv[0]*f_second_deriv.data()[3];

    }

    template<typename T>
    bool compute_direction(VectorX<T>& p, const Block<T>* b,
    const Block<T>* b2, VectorX<T>& first_deriv, T& gradient_magnitude, T& second_deriv_tangent)
    {
        if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
        {
            return false;
        }

        compute_gradient_h(p,b,b2,first_deriv,gradient_magnitude,second_deriv_tangent,true);
        first_deriv.normalize();
        return true;
    }


    template<typename T>    //normalized first derivatie
    bool compute_direction(VectorX<T>& p, const Block<T>* b, const Block<T>* b2, VectorX<T>& first_deriv)
    {
        if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
        {
            return false;
        }
        T gradient_magnitude; T second_deriv_tangent;
        compute_gradient_h(p,b,b2, first_deriv,gradient_magnitude,second_deriv_tangent);


        first_deriv.normalize();
        return true;
    }


    template<typename T>
    bool RKF45(T h, const Block<T>* b, 
    const Block<T>* b2, VectorX<T>& inital, VectorX<T>& result, T& gradient_magnitude, T& second_deriv_tangent, bool clock_wise)
    {
        // std::cout<<"ctrl_pts size "<<mfa_data->tmesh.tensor_prods[0].ctrl_pts.size()<<std::endl;
        
        VectorX<T> direction = VectorX<T>::Zero(inital.size());


        if(inital.size()==2)
        {
            if(clock_wise)
            {
                VectorX<T> v0 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(inital,b, b2,direction,gradient_magnitude,second_deriv_tangent))
                {
                    return false;
                }
                v0[0]=direction[1];v0[1]=-direction[0];

                VectorX<T> Y = inital + 0.5*h*v0;
                VectorX<T> v1 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(Y,b, b2,direction))
                {
                    return false;
                }
                v1[0]=direction[1];v1[1]=-direction[0];
                Y = inital + 0.5*h*v1;
                VectorX<T> v2 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(Y,b,b2,direction))
                {
                    return false;
                }
                v2[0]=direction[1];v2[1]=-direction[0];
                Y = inital + h*v2;
                VectorX<T> v3 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(Y,b, b2,direction))
                {
                    return false;
                }
                v3[0]=direction[1];v3[1]=-direction[0];
                result = inital + h/6*(v0+2*v1+2*v2+v3);

            }
            else
            {
                VectorX<T> v0 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(inital,b,b2,direction,gradient_magnitude,second_deriv_tangent))
                {
                    return false;
                }
                v0[0]=-direction[1];v0[1]=direction[0];
                VectorX<T> Y = inital + 0.5*h*v0;
                VectorX<T> v1 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(Y,b, b2,direction))
                {
                    return false;
                }
                v1[0]=-direction[1];v1[1]=direction[0];
                Y = inital + 0.5*h*v1;
                VectorX<T> v2 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(Y,b, b2,direction))
                {
                    return false;
                }
                v2[0]=-direction[1];v2[1]=direction[0];
                Y = inital + h*v2;
                VectorX<T> v3 = VectorX<T>::Zero(inital.size());
                if(!compute_direction(Y,b, b2,direction))
                {
                    return false;
                }
                v3[0]=-direction[1];v3[1]=direction[0];
                result = inital + h/6*(v0+2*v1+2*v2+v3);
                // std::cout<<"result "<<inital.transpose()<< result.transpose()<<" "<<v0.transpose()<<" "<<v1.transpose()<<" "<<v2.transpose()<<" "<<v3.transpose()<<std::endl;
            }


            if (!utility::In_Domain(result,b->core_mins,b->core_maxs))
            {
                return false;
            }
        }

        return true;
    } 


  

    template<typename T>
    // to save computation, compute gradient and hessian at the same time
    bool gradient_descent_compute_h_gradient_h(VectorX<T>& p, const Block<T>* b,
    const Block<T>* b2, T& h, VectorX<T>& h_first_deriv,T value, T threshold)
    {
        VectorX<T> f_first_deriv(p.size());
        compute_gradient(p,b,f_first_deriv);
        VectorX<T> g_first_deriv(p.size());
        compute_gradient(p,b2,g_first_deriv);
        h = f_first_deriv[0]*g_first_deriv[1]-f_first_deriv[1]*g_first_deriv[0];

        if(abs(h-value)<threshold)
        {
            return true;
        }

        MatrixX<T> f_second_deriv(p.size(),p.size());
        compute_hessian(p,b,f_second_deriv);
        
        MatrixX<T> g_second_deriv(p.size(),p.size());
        compute_hessian(p,b2,g_second_deriv);

        h_first_deriv[0]=f_first_deriv[0]*g_second_deriv.data()[1]+g_first_deriv[1]*f_second_deriv.data()[0]-f_first_deriv[1]*g_second_deriv.data()[0]-g_first_deriv[0]*f_second_deriv.data()[1];

        h_first_deriv[1]=f_first_deriv[0]*g_second_deriv.data()[3]+g_first_deriv[1]*f_second_deriv.data()[1]-f_first_deriv[1]*g_second_deriv.data()[1]-g_first_deriv[0]*f_second_deriv.data()[3];

        return false;


    }

    template<typename T>
    //Here value should input target function value
    bool gradient_descent(VectorX<T> &initial, const Block<T>* b, const Block<T>* b2, T value, int max_itr, T threshold)
    {
        VectorX<T> first_deriv(initial.size());


        T func_value;
        int i=0;


        while (i<max_itr)
        {

            if(gradient_descent_compute_h_gradient_h(initial,b,b2,func_value,first_deriv,value,threshold))
            {
                return true;
            }

            // int j=0;
            // while((initial-previous_point).squaredNorm()<4e-4*step_size_square)
            // {
            //     initial += direction;
            //     if(!utility::InDomain(initial))
            //     {
            //         return false;
            //     }
            //     // move out based on the previous direction
            //     compute_h_gradient_h(initial,mfa,mfa_data,func_value,first_deriv);
            // }
         

            T step_size = func_value/first_deriv.squaredNorm();


            initial -= step_size*first_deriv;

            if(!utility::In_Domain(initial,b->core_mins,b->core_maxs))
            {
                return false;
            }

            i++;       

        }
        return false;
       
    }

    template<typename T>
    bool correct_result(VectorX<T> &initial, const Block<T>* b, const Block<T>* b2,T ori_value, T threshold, int correction_max_itr)
    {
        T value;

        // if(!utility::InDomain(initial,domain_min,domain_range))
        // {
        //     return false;   
        // }

        compute_h(initial,b,b2,value);

        if(abs(value-ori_value)<threshold)
        {
            return true;
        }

        // VectorX<T> direction;
        // if(!record_result.empty())
        // {
        //     direction =0.25*(initial - record_result.back());
        // }
        // else
        // {
        //     direction = span_center-initial;
        // }

        return gradient_descent(initial,b,b2,ori_value,correction_max_itr,threshold);


    }


    // true means a trace, false means a loop
    template<typename T>
    bool tracing_one_direction_isocontour(T h,int max_step, const Block<T>* b,
    const Block<T>* b2, VectorX<T>& initial, std::vector<VectorX<T>>& result, std::vector<T>& gradient_magnitude, std::vector<T>& second_deriv_tangent, bool clock_wise, std::vector<std::vector<T>>& span_range, T ori_func_value, T threshold_correction, int correction_max_itr, T squared_threshold_distance_stop_tracing) //clock_wise: the direction to trace, span_range: [[min0,max0],[min1,max1]]
    {
        result.clear();
        gradient_magnitude.clear();

        VectorX<T> p_new = VectorX<T>::Zero(initial.size());
        VectorX<T> p_old = initial;

        T temp_gradient_magnitude;
        T temp_second_deriv_tangent;

        bool valid = RKF45(h,b,b2,initial,p_new,temp_gradient_magnitude,temp_second_deriv_tangent, clock_wise);

        result.emplace_back(initial);
        // gradient_magnitude.emplace_back(temp_gradient_magnitude);
        // second_deriv_tangent.emplace_back(temp_second_deriv_tangent);

        bool work_on_tracing = true;

        if(!valid)
        {
            // result.emplace_back(initial);
            // gradient_magnitude.emplace_back(temp_gradient_magnitude);
            // second_deriv_tangent.emplace_back(temp_second_deriv_tangent);
            work_on_tracing = false;
        }

    

        int step=1;


        // VectorX<T> span_center = VectorX<T>::Zero(initial.size());
        // for(int i=0;i<span_center.size();i++)
        // {
        //     span_center[i] = 0.5*(span_range[i][0]+span_range[i][1]);
        // }



        if(valid)
        {
            if(!correct_result(p_new,b,b2,ori_func_value,threshold_correction,correction_max_itr))
            {
                work_on_tracing = false;
            }
        }



        if(work_on_tracing)
        {
            T standard=(p_new-initial).squaredNorm();        
            if(standard<squared_threshold_distance_stop_tracing*h*h)
            {
                result.emplace_back(p_new);
                return true;
            }

            standard *= 1.21;


            while(utility::InBlock(span_range,p_new))
            {            
                if(((p_new-initial).squaredNorm()<standard || step>max_step) && step>1)
                {
                     result.emplace_back(p_new);
                    if(step>2)
                    {
                        return false;
                    }
                    else
                    {
                        return true;
                    }
                }
                result.emplace_back(p_new);
                p_old = p_new;

                // std::cout<<"step "<<step<<" "<<p_new.transpose()<<std::endl;

                valid = RKF45(h,b,b2,p_old,p_new,temp_gradient_magnitude,temp_second_deriv_tangent,clock_wise);
                // gradient_magnitude.emplace_back(temp_gradient_magnitude);
                // second_deriv_tangent.emplace_back(temp_second_deriv_tangent);


                if(!valid)
                {
                    break;
                }  

                if(!correct_result(p_new,b,b2,ori_func_value,threshold_correction,correction_max_itr))
                {
                    break;
                }


                if((p_new-p_old).squaredNorm()<squared_threshold_distance_stop_tracing*h*h)
                {
                     result.emplace_back(p_new);
                    return true;
                }

                step++; 
                
            }

        }


        // std::cout<<"finish step "<<valid<<std::endl;
        // std::cout<<p_new.transpose()<<std::endl;


        // use RKF45 one more time, with half step size
        valid = RKF45(0.5*h,b,b2,p_old,p_new,temp_gradient_magnitude,temp_second_deriv_tangent,clock_wise);



        if(valid)
        {
            if(correct_result(p_new,b,b2,ori_func_value,threshold_correction,correction_max_itr))
            {
                if(utility::InBlock(span_range,p_new))
                {
                    result.emplace_back(p_new);
                    // VectorX<T> gradient_ = p_new;
                    // compute_gradient_h(p_new,mfa,mfa_data,mfa2,mfa_data2,gradient_,temp_gradient_magnitude,temp_second_deriv_tangent,domain_min,domain_range,true);
                    // gradient_magnitude.emplace_back(temp_gradient_magnitude);
                    // second_deriv_tangent.emplace_back(temp_second_deriv_tangent);
                    p_old = p_new;
                }
            }
        }

        valid = RKF45(0.1*h,b,b2,p_old,p_new,temp_gradient_magnitude,temp_second_deriv_tangent,clock_wise);

        if(valid)
        {
            if(correct_result(p_new,b,b2,ori_func_value,threshold_correction,correction_max_itr))
            {
                if(utility::InBlock(span_range,p_new))
                {                
                    result.emplace_back(p_new);
                    
                }
            }
        }

        return true;
        
    }

    
    template<typename T>
    bool tracing_isocontour(T h,int max_step, const Block<T>* b, const Block<T>* b2, VectorX<T>& initial, std::vector<VectorX<T>>& result, std::vector<T>& gradient_magnitude,
    std::vector<T>& second_deriv_tangent, std::vector<std::vector<T>>& span_range, T ori_func_value, T threshold_correction, int correction_max_itr,T squared_threshold_distance_stop_tracing) //clock_wise: the direction to trace, span_range: [[min0,max0],[min1,max1]], return false for a loop
    {
        result.clear();
        gradient_magnitude.clear();
        second_deriv_tangent.clear();


        if(tracing_one_direction_isocontour(h,max_step,b,b2,initial,result,gradient_magnitude,second_deriv_tangent,true,span_range,ori_func_value,threshold_correction,correction_max_itr,squared_threshold_distance_stop_tracing))
        {
            std::vector<VectorX<T>> temp_result;
            std::vector<T> temp_gradient_magnitude;        
            std::vector<T> temp_second_deriv_tangent;   

            tracing_one_direction_isocontour(h,max_step,b,b2, initial,temp_result,temp_gradient_magnitude,temp_second_deriv_tangent,false,span_range,ori_func_value,threshold_correction,correction_max_itr,squared_threshold_distance_stop_tracing);
            std::reverse(result.begin(),result.end());
            // std::reverse(gradient_magnitude.begin(),gradient_magnitude.end());


            // for(int i=0;i<result.size();i++)
            // {
            //     std::cout<<result[i].transpose()<<" "<<gradient_magnitude[i]<<std::endl;
            // }
            
            // std::cout<<"-----------------"<<std::endl;

            // for(int i=0;i<temp_result.size();i++)
            // {
            //     std::cout<<temp_result[i].transpose()<<" "<<temp_gradient_magnitude[i]<<std::endl;
            // }


            if(temp_result.size()>1)
            {
                result.pop_back();
                // gradient_magnitude.pop_back();
                // second_deriv_tangent.pop_back();

                // for(int i=0;i<temp_result.size();i++)
                // {
                //     std::cout<<temp_result[i].transpose()<<" "<<temp_gradient_magnitude[i]<<std::endl;
                // }


                result.insert(result.end(),temp_result.begin(),temp_result.end());
                // gradient_magnitude.insert(gradient_magnitude.end(),temp_gradient_magnitude.begin(),temp_gradient_magnitude.end());
                // second_deriv_tangent.insert(second_deriv_tangent.end(),temp_second_deriv_tangent.begin(),temp_second_deriv_tangent.end());
                
            }
            return true;
        }        

        return false;
    }

    template<typename T>
    int index_map(int length_A, int length_B, int index_A)
    {
        return (int)((T)(length_B)/(T)(length_A)*(T)index_A);
    }


    template<typename T>
    //check if two traces are the same. The step sizes of two traces are the same
    bool trace_duplicate(std::vector<Eigen::VectorX<T>>& trace1,std::vector<Eigen::VectorX<T>>& trace2, T threshold_square)
    {
        int corresponding_index=0;
        T distance_2=(trace1[0]-trace2[0]).squaredNorm();
        T temp_dis;
        for (int i = 1; i < trace2.size(); i++)
        {
            temp_dis = (trace1[0]-trace2[i]).squaredNorm();
            if (temp_dis < distance_2)
            {
                distance_2 = temp_dis;
                corresponding_index = i;
            }
        }

        if(distance_2 > threshold_square)
        {
            return false;
        }

        int mid_corresponding_index = (corresponding_index+(trace2.size()>>1))%trace2.size();
        int mid_1_index=trace1.size()>>1;

        distance_2 = (trace1[mid_1_index]-trace2[mid_corresponding_index]).squaredNorm();
        for(int i=-1;i<2;i+=2)
        {
            temp_dis = (trace1[mid_1_index]-trace2[(mid_corresponding_index+i)%trace2.size()]).squaredNorm();
            if(temp_dis<distance_2)
            {
                distance_2 = temp_dis;
            }
        }
        if(distance_2 > threshold_square)
        {
            return false;
        }


        return true;

    }


    // template<typename T>
    // //check if two traces are the same. The step sizes of two traces are the same
    // bool trace_duplicate(std::vector<Eigen::VectorX<T>>& trace1,std::vector<Eigen::VectorX<T>>& trace2, T threshold_square)
    // {
    //     int corresponding_index=0;
    //     T distance_2=(trace1[0]-trace2[0]).squaredNorm();
    //     T temp_dis;
    //     for (int i = 1; i < trace2.size(); i++)
    //     {
    //         temp_dis = (trace1[0]-trace2[i]).squaredNorm();
    //         if (temp_dis < distance_2)
    //         {
    //             distance_2 = temp_dis;
    //             corresponding_index = i;
    //         }
    //     }

    //     if(distance_2 > threshold_square)
    //     {
    //         return false;
    //     }
    //     // point 0 in trace1 is corresponding to point corresponding_index in trace2
    //     // the step size are fixed, traces 1 and 2 has the same step size
    //     // But it's possible that trace 1 and 2 may have different steps. 
    //     // Find last point in trace1 that is corresponding to a point in trace2 based on the relation
    //     //
    //     temp_dis = 0;
    //     T last_distance_2 = (trace1.back()-trace2[corresponding_index]).squaredNorm();
    //     int last_point_corresponding_index = corresponding_index;
    //     temp_dis = (trace1.back()-trace2[(corresponding_index-1)%trace2.size()]).squaredNorm();
    //     if (temp_dis < last_distance_2)
    //     {
    //         last_distance_2 = temp_dis;
    //         last_point_corresponding_index = (corresponding_index-1)%trace2.size();
    //     }
        
    //     if(last_distance_2>distance_2)
    //     {
    //         distance_2 = last_distance_2;
    //         if(last_distance_2 > threshold_square)
    //         {
    //             return false;
    //         }
    //     }
        
    //     int record_corresponding_index;
    //     for (int i=1;i<trace1.size()-1;i++)
    //     {
    //         record_corresponding_index = index_map<T>(trace1.size()-1,trace2.size()-(int)(last_point_corresponding_index!=corresponding_index),i);
    //         temp_dis = (trace1[i]-trace2[record_corresponding_index]).squaredNorm();
    //         T temp_dis_2 = (trace1[i]-trace2[(record_corresponding_index-1)%trace2.size()]).squaredNorm();
    //         if (temp_dis > temp_dis_2)
    //         {
    //             temp_dis = temp_dis_2;
    //         }
    //         temp_dis_2 = (trace1[i]-trace2[(record_corresponding_index+1)%trace2.size()]).squaredNorm();
    //         if (temp_dis > temp_dis_2)
    //         {
    //             temp_dis = temp_dis_2;
    //         }
    //         if(temp_dis>distance_2)
    //         {
    //             distance_2 = temp_dis;
    //             if(temp_dis > threshold_square)
    //             {
    //                 return false;
    //             }
    //         }

    //     }

    //     return true;


    // }
    



//use this threshold to locate some areas that contains zero magnitude, then chose a point from it
    template<typename T>
    void findZeroArea(std::vector<T>&gradient_magnitude, T threshold_square, std::vector<int>& split_index)
    {
        bool start_recording_zero=false;
        T record_min_magitude = std::numeric_limits<T>::max();
        int record_min_index = -1;


        int start_index=0;
        if(gradient_magnitude[0]<threshold_square)
        {
            for(int i=1;i<gradient_magnitude.size();++i)
            {
                if(gradient_magnitude[i]>threshold_square)
                {
                    start_index = i;
                    break;
                }
            }

            if(start_index==0)
            {
                return;
            }
        }

        int end_index=gradient_magnitude.size()-1;
        if(gradient_magnitude[end_index]<threshold_square)
        {
            for (int i = end_index-1; i > start_index; i--)
            {
                if(gradient_magnitude[i]>threshold_square)
                {
                    end_index = i;
                    break;
                }
            }

            if(end_index==gradient_magnitude.size()-1)
            {
                return;
            }
        }


        for(int i=start_index;i<=end_index;++i)
        {
            if(!start_recording_zero)
            {
                if(gradient_magnitude[i]<threshold_square)
                {
                    start_recording_zero=true;
                    record_min_magitude = std::numeric_limits<T>::max();
                }
            }
            else
            {
                if(gradient_magnitude[i]>threshold_square)
                {
                    start_recording_zero=false;
                    split_index.emplace_back(record_min_index);
                }
            
            }

            if(start_recording_zero)
            {
                if(record_min_magitude>gradient_magnitude[i])
                {
                    record_min_magitude = gradient_magnitude[i];
                    record_min_index = i;
                }
            }

        }

    }

    template<typename T>
    void split_trace(std::vector<VectorX<T>>& trace, std::vector<T>&gradient_magnitude,
    std::vector<T>& second_deriv_tangent,
    std::vector<std::vector<VectorX<T>>>& result,
    std::vector<std::vector<T>>& gradient_magnitude_result,
    std::vector<std::vector<T>>& second_deriv_tangent_result,
    T threshold_square)
    {
        std::vector<int> split_index;

        // findZeroArea(gradient_magnitude,threshold_square,split_index);


        // std::cout<<threshold_square<<std::endl;
        // std::cout<<"split index "<<std::endl;
        // for(int i=0;i<gradient_magnitude.size();++i)
        // {
        //     std::cout<<gradient_magnitude[i]<<" ";
        // }
        // std::cout<<std::endl;
        // for(int i=0;i<split_index.size();++i)
        // {
        //     std::cout<<split_index[i]<<" ";
        // }
        // std::cout<<std::endl;

        result.clear();
        // gradient_magnitude_result.clear();
        // second_deriv_tangent_result.clear();

        // if(split_index.empty())
        // {
            result.emplace_back(trace);
            // gradient_magnitude_result.emplace_back(gradient_magnitude);
            // second_deriv_tangent_result.emplace_back(second_deriv_tangent);
            // return;
        // }

        ////the vector [a-1,b] [b-1,c] ...
        // result.emplace_back(trace.begin(),trace.begin()+split_index[0]);
        // gradient_magnitude_result.emplace_back(gradient_magnitude.begin(),gradient_magnitude.begin()+split_index[0]);
        // second_deriv_tangent_result.emplace_back(second_deriv_tangent.begin(),second_deriv_tangent.begin()+split_index[0]);
        
        // for(int i=1;i<split_index.size();i++)
        // {
        //     result.emplace_back(trace.begin()+split_index[i-1],trace.begin()+split_index[i]);
        //     gradient_magnitude_result.emplace_back(gradient_magnitude.begin()+split_index[i-1],gradient_magnitude.begin()+split_index[i]);
        //     second_deriv_tangent_result.emplace_back(second_deriv_tangent.begin()+split_index[i-1],second_deriv_tangent.begin()+split_index[i]);
        // }

        // result.emplace_back(trace.begin()+split_index.back(),trace.end());
        // gradient_magnitude_result.emplace_back(gradient_magnitude.begin()+split_index.back(),gradient_magnitude.end());
        // second_deriv_tangent_result.emplace_back(second_deriv_tangent.begin()+split_index.back(),second_deriv_tangent.end());
        
    }

        //check if split index already exists
    bool check_index_already_exist(std::vector<int>& split_index, int current_split_index)
    {
        if(split_index.empty())
        {
            split_index.emplace_back(current_split_index);
            return true;
        }
        else
        {
            bool all_not_unique = true;
            for(int k=0;k<split_index.size();++k)
            {
                if(std::abs(split_index[k] -current_split_index)<2)
                {
                    all_not_unique = false;
                    break;
                }
            }
            if(all_not_unique)
            {
                split_index.emplace_back(current_split_index);
                return true;
            }          
        }       
        return false;
    }


    //check if two traces share the same part then go to different directions
    template<typename T>
    bool check_index_share_same_part(std::vector<int>& split_index_i, std::vector<int>& split_index_j, std::vector<VectorX<T>>& vec_i,std::vector<VectorX<T>>& vec_j, bool same_start,T threshold_square)
    {
        if(vec_i.size()<3 || vec_j.size()<3)
        {
            return false;
        }
        int check_index_i = 1;
        int check_index_j = 1;



        if(!same_start)
        {
            check_index_i = vec_i.size()-2;
            check_index_j = vec_j.size()-2;

        }

        T k0 = (vec_i[check_index_i]-vec_j[check_index_j]).squaredNorm();
        if(k0>threshold_square)
        {
            return false;
        }

        int update_unit=1;
        int loop_end_i= vec_i.size();
        int loop_end_j= vec_j.size();
        

        if(!same_start)
        {
            update_unit = -1;
            loop_end_i = -1;
            loop_end_j = -1;
        }


        int split_index_i_new=-1;
        int split_index_j_new=-1;
        if(same_start)
        {            
            for(int i=check_index_i+1,j=check_index_j+1;i<loop_end_i,j<loop_end_j;i++,j++)
            {            
                if((vec_i[i]-vec_j[j]).squaredNorm()>threshold_square)
                {
                    split_index_i_new = i;
                    split_index_j_new = j;
                    break;                    
                }
            }
        }
        else
        { 
            for(int i=check_index_i-1,j=check_index_j-1;i>loop_end_i,j>loop_end_j;i--,j--)
            {            
                if((vec_i[i]-vec_j[j]).squaredNorm()>threshold_square)
                {
                    split_index_i_new = i+1;
                    split_index_j_new = j+1;
                    break;
                }
            }

        }

        if(split_index_i_new!=-1)
        { 
            check_index_already_exist(split_index_i,split_index_i_new);
            check_index_already_exist(split_index_j,split_index_j_new);
            return true;
        }
        return false;
    }


    template<typename T>
    bool check_vector_part_of_loop(int i, std::vector<VectorX<T>>& vector_loop, int j, std::vector<VectorX<T>>& vector_j, T threshold_square)
    {
        T k0 = (vector_loop[0]-vector_j[0]).squaredNorm();
        int index = -1;
        for(int k=0;k<vector_loop.size();++k)
        {
            if((vector_loop[k]-vector_j[0]).squaredNorm()<threshold_square)
            {
                index = k;
                break;
            }
        }

        int index2 = -1;

        for(int k=0;k<vector_loop.size();++k)
        {
            if((vector_loop[k]-vector_j.back()).squaredNorm()<threshold_square)
            {
                index2 = k;
                break;
            }
        }


        if(index2>-1 && index>-1)
        {
            return true;
        }
    
        return false;

    }


    // if a vector is a part of another vector, then keep the short vector and divide the long vector into two vectors based on the short vector and remove the part that is the same as the short vector
    //remove_info is consistent with split_index, to show we should remove the front or back part of the vector from split_index. For the case that the entire vector should be removed, the corresponding element in remove_info is -1

    template<typename T>
    bool check_part_of_vector(int i, std::vector<VectorX<T>>& vector_i, int j, std::vector<VectorX<T>>& vector_j, T threshold_square, std::vector<std::vector<int>>& split_index)
    {

        T k0 = (vector_i[0]-vector_j[0]).squaredNorm();
        T k1 = (vector_i.back()-vector_j.back()).squaredNorm();
        if(k0>threshold_square && k1>threshold_square)
        {
            return false;
        }
        else if(k0<threshold_square)
        {
            bool i_short=true;
            T k2;
            if(vector_i.size()>vector_j.size())
            {
                i_short=false;
                k2= (vector_i[vector_j.size()-1]-vector_j.back()).squaredNorm();        
            }
            else
            {
                k2= (vector_j[vector_i.size()-1]-vector_i.back()).squaredNorm();
            }
            if(k2>threshold_square)
            {
                return check_index_share_same_part(split_index[i],split_index[j],vector_i,vector_j,true,threshold_square);
            }
            if(i_short)
            {
                return check_index_already_exist(split_index[j],vector_i.size());
            }
            else
            {
                return check_index_already_exist(split_index[i],vector_j.size());
            }
        }
        else//(k1<threshold_square)
        {
            bool i_short=true;
            T k2;
            if(vector_i.size()>vector_j.size())
            {
                i_short=false;
                k2= (vector_i[vector_i.size() - vector_j.size()]-vector_j[0]).squaredNorm();        
            }
            else
            {
                k2= (vector_j[vector_j.size()-vector_i.size()]-vector_i[0]).squaredNorm();
            }
            if(k2>threshold_square)
            {
                return check_index_share_same_part(split_index[i],split_index[j],vector_i,vector_j,false,threshold_square);
            }
            if(i_short)
            {
                return check_index_already_exist(split_index[j],(int)(vector_j.size())-(int)(vector_i.size()));
            }
            else
            {
                return check_index_already_exist(split_index[i],(int)(vector_i.size())-(int)(vector_j.size()));
            }
        }

        return false;
    }


    template<typename T> 
    void split_vec_on_index(std::vector<int>& split_index, std::vector<std::vector<VectorX<T>>>& new_trace, std::vector<bool>& new_loop,std::vector<VectorX<T>>& ori_trace, bool ori_loop)
    {
        new_trace.emplace_back(ori_trace.begin(),ori_trace.begin()+split_index[0]);
        new_loop.emplace_back(ori_loop);
        for(int i=1;i<split_index.size();++i)
        {
            new_trace.emplace_back(ori_trace.begin()+split_index[i-1],ori_trace.begin()+split_index[i]);
            new_loop.emplace_back(ori_loop);
        }
        new_trace.emplace_back(ori_trace.begin()+split_index.back(),ori_trace.end());
        new_loop.emplace_back(ori_loop);
    }




    //deduplicate
    template<typename T>
    void deduplicate(std::vector<std::vector<VectorX<T>>>& record_trace, std::vector<bool>& record_loop, std::vector<std::vector<VectorX<T>>>& trace, std::vector<bool>& loop, T threshold_square)
    {        
        record_trace.clear();
        record_loop.clear();

        for(int i=0;i<trace.size();++i)
        {
            bool find_same = false;
            if(!loop[i])
            {
                if(trace[i].size()>1)
                {
                    for(int j=0;j<record_trace.size();++j)
                    {
                        if(!record_loop[j])
                        {
                            if((trace[i][0]-record_trace[j][0]).squaredNorm()<threshold_square && (trace[i].back()-record_trace[j].back()).squaredNorm()<threshold_square)
                            {
                                find_same = true;
                                break;
                            }
                        }
                    }
                }
                else if(trace[i].empty())
                {
                    find_same = true;
                }
            }
            if(!find_same)
            {
                record_trace.emplace_back(trace[i]);
                record_loop.emplace_back(loop[i]);
            }
        }
    }

    template<typename T>
    void split_and_deduplicate_based_on_index(TraceInSpan<T>& trace_in_span, std::vector<std::vector<int>>& split_index,T threshold_square)
    {
        TraceInSpan<T> temp_trace_in_span;
        temp_trace_in_span.traces.reserve(trace_in_span.traces.size());
        temp_trace_in_span.is_loop.reserve(trace_in_span.traces.size());


        for(auto i=0;i<trace_in_span.traces.size();++i)
        {
            if(!trace_in_span.traces[i].empty())
            {
                if(split_index[i].empty())
                {
                    temp_trace_in_span.traces.emplace_back(trace_in_span.traces[i]);
                    temp_trace_in_span.is_loop.emplace_back(trace_in_span.is_loop[i]);
                }
                else
                {
                    split_vec_on_index(split_index[i],temp_trace_in_span.traces,temp_trace_in_span.is_loop,trace_in_span.traces[i],trace_in_span.is_loop[i]);
                }
            }
        }

        // if(!temp_trace_in_span.traces.empty())
        // {
        //     std::cout<<temp_trace_in_span.traces[0].size()<<" "<<temp_trace_in_span.traces.size()<<std::endl;
        // }

        //deduplicate
        deduplicate(trace_in_span.traces,trace_in_span.is_loop,temp_trace_in_span.traces,temp_trace_in_span.is_loop,threshold_square);
        
    }


    template<typename T>
    bool check_a_point_is_in_trace(VectorX<T>& point, std::vector<VectorX<T>>& vector_j, T threshold_square)
    {
        for(auto i=0;i<vector_j.size();++i)
        {
            if((vector_j[i]-point).squaredNorm()<threshold_square)
            {
                if((vector_j[(i+1)%vector_j.size()]-point).squaredNorm()<threshold_square)
                {
                    return true;
                }
            }
        }
        return false;
    }


    //check trace is a part of a loop
    template<typename T>
    void check_duplication_3(TraceInSpan<T>& trace_in_span,T threshold_square)
    {


        std::vector<bool> remove_trace_index(trace_in_span.traces.size(),false);

            
        for(auto i=0;i<trace_in_span.traces.size();++i)
        {
            if(trace_in_span.is_loop[i])
            {
                for(auto j=0;j<trace_in_span.traces.size();++j)
                {
                    if((!trace_in_span.is_loop[j]) &&(!remove_trace_index[j]))
                    {
                        if(check_vector_part_of_loop(i,trace_in_span.traces[i],j,trace_in_span.traces[j],threshold_square))
                        {
                            remove_trace_index[j]=true;
                        }
                    }
                }
            }
        }


        for(auto i=0;i<trace_in_span.traces.size();++i)
        {
            if(trace_in_span.traces[i].size()==1 && (!remove_trace_index[i]))
            {
                for(auto j=0;j<trace_in_span.traces.size();++j)
                {
                    if(trace_in_span.traces[j].size()>1 && (!remove_trace_index[j]))
                    {
                        if(check_a_point_is_in_trace(trace_in_span.traces[i][0],trace_in_span.traces[j],threshold_square))
                        {
                            remove_trace_index[i]=true;
                            break;
                        }                        
                    }
                }
                
            }

        }


        bool all_false = std::all_of(remove_trace_index.begin(), remove_trace_index.end(), [](bool x) { return !x; });

        if(all_false)
        {
            return;
        }

        TraceInSpan<T> temp_trace_in_span;
        temp_trace_in_span.traces.reserve(trace_in_span.traces.size());
        temp_trace_in_span.is_loop.reserve(trace_in_span.traces.size());
        for(auto i=0;i<trace_in_span.traces.size();++i)
        {
            if(!remove_trace_index[i])
            {
                temp_trace_in_span.traces.emplace_back(trace_in_span.traces[i]);
                temp_trace_in_span.is_loop.emplace_back(trace_in_span.is_loop[i]);
            }
        }

        trace_in_span = temp_trace_in_span;
        
    }

    // check is a trace is a part of other trances
    template<typename T>
    void check_duplication_2(TraceInSpan<T>& trace_in_span,T threshold_square)
    {
        std::vector<std::vector<int>> split_index(trace_in_span.traces.size());
            
        for(auto i=0;i<trace_in_span.traces.size();++i)
        {
            if(!trace_in_span.is_loop[i])
            {
                for(auto j=i+1;j<trace_in_span.traces.size();++j)
                {
                    if(!trace_in_span.is_loop[j])
                    {
                        check_part_of_vector(i,trace_in_span.traces[i],j,trace_in_span.traces[j],threshold_square,split_index);
                    }
                }
            }
        }


        bool all_empty = std::all_of(split_index.begin(), split_index.end(), [](std::vector<int>& x) { return x.empty(); });

        if(all_empty)
        {
            return;
        }

        for(auto i=split_index.begin();i!=split_index.end();++i)
        {
            if(!i->empty())
            {
                std::sort(i->begin(),i->end());
            }
        }

        // // //test
        // for(int i=0;i<split_index.size();++i)
        // {
        //     if(split_index[i].empty())
        //     {
        //         continue;
        //     }
        //     for(int j=1;j<split_index[i].size();++j)
        //     {
        //         if(split_index[i][j]-split_index[i][j-1]<1)
        //         {
        //             std::cout<<"error "<<i<<" "<<j<<" "<<split_index[i][j]<<" "<<split_index[i][j-1]<<std::endl;
        //         }
        //     }
        //     if(split_index[i][0]==0)
        //     {
        //         std::cout<<"error START"<<i<<" "<<" "<<trace_in_span.traces[i].size()<<std::endl;
        //     }
        //     // if(split_index[i].back() == trace_in_span.traces[i].size()-1)
        //     // {
        //     //     std::cout<<"error END"<<i<<" "<<" "<<trace_in_span.traces[i].size()<<std::endl;
        //     // }

        // }


        split_and_deduplicate_based_on_index(trace_in_span,split_index,threshold_square);
        
    }


    template<typename T>
    bool check_duplication(TraceInSpan<T>& trace_in_span, std::vector<VectorX<T>>& result, bool is_loop, T threshold_square)
    {
        if(is_loop)
        {
            for(auto i=0;i<trace_in_span.traces.size();++i)
            {
                if(trace_in_span.is_loop[i])
                {
                    if(trace_duplicate(result,trace_in_span.traces[i],threshold_square))
                    {
                        return true;
                    }
                }
            }
        }
        else
        {           
                // std::cout<<"trace_threshold_square "<<threshold_square<<std::endl;
            // std::cout<<"check duplication "<<threshold_square<<std::endl;
            if(result.empty())
            {
                return true;
            }
            for(auto i=0;i<trace_in_span.traces.size();++i)
            {
                if(!trace_in_span.is_loop[i])
                {

                    if((result[0]-trace_in_span.traces[i][0]).squaredNorm()<threshold_square)
                    {
                        if((result.back()-trace_in_span.traces[i].back()).squaredNorm()<threshold_square)
                        {                            
                            return true;
                        }
                    }                    
                                       
                }
            }
            
        }
        return false;
    }


    template<typename T>
    void tracing_isocontour(T h,int max_step, const Block<T>* b, 
    const Block<T>* b2, std::vector<VectorX<T>>& initial, T trace_threshold_square,std::vector<std::vector<T>>& span_range,TraceInSpan<T>& trace_in_span,T ori_func_value, T threshold_correction, T trace_split_grad_square_threshold ,int correction_max_itr,T squared_threshold_distance_stop_tracing)
    {
        std::vector<VectorX<T>> result;
        std::vector<T> record_gradient_magnitude;
        std::vector<T> record_second_deriv_tangent;

        for (int i = 0; i < initial.size(); ++i)
        {

            bool not_loop= tracing_isocontour(h,max_step,b,b2,initial[i],result,record_gradient_magnitude,record_second_deriv_tangent,span_range,ori_func_value,threshold_correction,correction_max_itr,squared_threshold_distance_stop_tracing);

            std::vector<std::vector<VectorX<T>>> split_result;
            // std::vector<std::vector<T>> split_gradient_magnitude_result;
            // std::vector<std::vector<T>> split_second_deriv_tangent_result;

            if(not_loop)
            {            
                split_result.emplace_back(result);

                // split_trace(result,record_gradient_magnitude,record_second_deriv_tangent,split_result,split_gradient_magnitude_result,split_second_deriv_tangent_result,trace_split_grad_square_threshold);
            }
            // std::cout<<"tracing "<<i<<" finished "<<not_loop<<std::endl;
            // std::cout<<"test "<<trace_in_span.traces.size()<<std::endl;
            // std::cout<<trace_in_span.traces.size()<<std::endl;
            // center = remove_duplicate_trace::compute_center(result);

            // for(int k=0;k<record_gradient_magnitude.size;++k)
            // {
            //     std::cout<<record_gradient_magnitude[k]<<" ";
            // }
            // std::cout<<std::endl;
            
            if(!not_loop)
            {
                if(!check_duplication(trace_in_span,result,!not_loop,trace_threshold_square))
                {
                    trace_in_span.traces.emplace_back(result);
                    // trace_in_span.gradient_magnitude.emplace_back(record_gradient_magnitude);
                    trace_in_span.is_loop.emplace_back(!not_loop);
                    // trace_in_span.second_deriv_tangent.emplace_back(record_second_deriv_tangent);


                }                
            }
            else
            {
                for(auto j=0;j<split_result.size();++j)
                {
                    // std::cout<<"split vector "<<std::endl;
                    // for(int k=0;k<split_gradient_magnitude_result[j].size();++k)
                    // {
                    //     std::cout<<split_gradient_magnitude_result[j][k]<<" ";
                    // }
                    // std::cout<<std::endl;
                    // if(split_result[j].size()<20)
                    if(!check_duplication(trace_in_span,split_result[j],!not_loop,trace_threshold_square))
                    {
                        trace_in_span.traces.emplace_back(split_result[j]);
                        // trace_in_span.gradient_magnitude.emplace_back(split_gradient_magnitude_result[j]);
                        trace_in_span.is_loop.emplace_back(!not_loop);
                        // trace_in_span.second_deriv_tangent.emplace_back(split_second_deriv_tangent_result[j]);


                    }                    
                    
                }
            }

        }

        check_duplication_2(trace_in_span,trace_threshold_square);
        check_duplication_2(trace_in_span,trace_threshold_square);
        check_duplication_3(trace_in_span,trace_threshold_square);
    }



    
    template<typename T>
    void test_points_feature_(const Block<T>* b, 
    const Block<T>* b2,TraceInSpan<T>& traces_in_span, std::vector<VectorX<T>>& f_gradient_, std::vector<VectorX<T>>& h_gradient_, 
    std::vector<T>& h_gradient_magnitude_, std::vector<T>& second_derive_tangent_, std::vector<VectorX<T>>& g_gradient_, std::vector<T>& h_)
    {
        VectorX<T> f_gradient(2);
        T f;
        VectorX<T> h_gradient(2);
        T h_gradient_magnitude;
        T second_derive_tangent;
        VectorX<T> g_gradient(2);
        MatrixX<T> hessian(2,2);
        T h;

        for(int i=0;i<traces_in_span.traces.size();++i)
        {
            for( int j=0;j<traces_in_span.traces[i].size();++j)
            {
                
                compute_gradient(traces_in_span.traces[i][j],b,f_gradient);
                compute_gradient_h(traces_in_span.traces[i][j],b,b2, h_gradient,h_gradient_magnitude,second_derive_tangent,true);
                compute_gradient(traces_in_span.traces[i][j],b2,g_gradient);
                compute_h(traces_in_span.traces[i][j],b,b2,h);
                f_gradient_.emplace_back(f_gradient);
                h_gradient_.emplace_back(h_gradient);
                g_gradient_.emplace_back(g_gradient);
                h_.emplace_back(h);
                h_gradient_magnitude_.emplace_back(h_gradient_magnitude);
                second_derive_tangent_.emplace_back(second_derive_tangent);
                
            }
            

        }
    }

    template<typename T>
    void test_points_feature(const Block<T>* b, const Block<T>* b2, std::vector<TraceInSpan<T>>& traces_in_span)
    {
        std::vector<VectorX<T>> f_gradient_;
        std::vector<VectorX<T>> h_gradient_;
        std::vector<T> f_gradient_magnitude_;
        std::vector<T> h_;
        std::vector<T> second_derive_tangent_;
        std::vector<VectorX<T>> g_gradient_;
    
        for(int i=0;i<traces_in_span.size();++i)
        {
            // if(abs(traces_in_span[i].traces[0][0][0])>0.875 && traces_in_span[i].traces[0][0][1]>0.875)
            // {
                test_points_feature_(b,b2,traces_in_span[i],f_gradient_,h_gradient_,f_gradient_magnitude_,second_derive_tangent_,g_gradient_,h_);
            // }
        }

        T max = *std::max_element(h_.begin(),h_.end());
        T min = *std::min_element(h_.begin(),h_.end());

        std::cout<<"max h "<<max<<std::endl;
        std::cout<<"min h "<<min<<std::endl;
        // for(int i=0;i<f_gradient_.size();++i)
        // {
        //     // std::cout<<"====="<<std::endl;
        //     std::cout<<h_[i]<<std::endl;
        //     // std::cout<<f_gradient_[i].transpose()<<std::endl;
        //     // std::cout<<g_gradient_[i].transpose()<<std::endl;
        //     // std::cout<<h_gradient_[i].transpose()<<std::endl;
        //     // std::cout<<f_gradient_magnitude_[i]<<std::endl;
        //     // std::cout<<second_derive_tangent_[i]<<std::endl;
        // }

    }


    template<typename T>
    void find_jacobi_set(T h,int max_step, const Block<T>* b, 
    const Block<T>* b2, std::vector<std::vector<VectorX<T>>>& initial,
    std::vector<size_t>& span_index, T trace_threshold_square, std::vector<TraceInSpan<T>>& traces_in_span, T threshold_correction, T trace_split_grad_square_threshold, int correction_max_itr,T squared_threshold_distance_stop_tracing)
    {
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);


        tbb::affinity_partitioner ap;
        auto domain_range = b->core_maxs-b->core_mins;

        // std::cout<<"range__ "<< domain_range.transpose()<<std::endl;
        // std::cout<<domain_min.transpose()<<std::endl;
        tbb::parallel_for(tbb::blocked_range<size_t>(0,initial.size()),[&](const tbb::blocked_range<size_t>& r)
        {
            for(auto i = r.begin(); i != r.end(); ++i)
            {        
                // std::cout<<"find_isocontour "<<initial[i].size()<<std::endl;
                VectorXi span_domain_index(b->mfa->var(0).p.size());
                //decode the span index in every dimension
                utility::obtainDomainIndex(span_index[i],span_domain_index,number_in_every_domain);
                // std::cout<<span_index[i]<<std::endl;
                // std::cout<<span_domain_index.transpose()<<std::endl;
                

                std::vector<std::vector<T>> span_range(number_in_every_domain.size());
                for(int j=0;j<span_domain_index.size();++j)
                {    
                    span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[j][span_domain_index[j]+b->mfa->var(0).p[j]]*domain_range[j]+b->core_mins[j]);
                    span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[j][span_domain_index[j]+1+b->mfa->var(0).p[j]]*domain_range[j]+b->core_mins[j]);
                }

                // if(span_range[0][0]<-2.14689 && span_range[0][1]>-2.14689 && 
                // span_range[1][0]<0.35746 && span_range[1][1]>0.35746)
                // {                
                //     std::cout<<span_range[0][0]<<std::endl;
                //     std::cout<<span_range[0][1]<<std::endl;
                //     std::cout<<span_range[1][0]<<std::endl;
                //     std::cout<<span_range[1][1]<<std::endl;
                    

                    tracing_isocontour(h,max_step,b,b2,initial[i],trace_threshold_square,span_range,traces_in_span[i],0.0,threshold_correction,trace_split_grad_square_threshold,correction_max_itr,squared_threshold_distance_stop_tracing);
                // }
            }
        },ap);           

    }



    // template<typename T>
    // void extract_ridge_valley_to_matrix(int split_ridge_valley, std::vector<TraceInSpan<T>>& traces_in_span,MatrixX<T>& ridge, MatrixX<T>& valley, MatrixX<T>& pseudo_ridge, MatrixX<T>& pseudo_valley, const VectorX<T>& core_min, const VectorX<T>& core_max)
    // {
    //     auto range = core_max-core_min;

    //     std::cout<<"test "<<std::endl;

    //     size_t point_size=0;
    //     for(auto i = 0;i<traces_in_span.size();++i)
    //     {
    //         point_size += traces_in_span[i].ridges.size();
    //     }
        

    //     // // std::cout<<point_size<<std::endl;

    //     // for(auto i = 0;i<traces_in_span.size();++i){

    //     //     std::cout<<traces_in_span[i].ridges.size()<<std::endl;
    //     //     std::cout<<traces_in_span[i].ridge_trace_index.size()<<std::endl;
        
    //     //     for(int j=0;j<traces_in_span[i].traces.size();j++)
    //     //     {
    //     //         std::cout<<"test "<<traces_in_span[i].traces[j].size()<<" "<< traces_in_span[i].traces_values[j].size()<<" ";
    //     //     }
    //     //     std::cout<<std::endl;
    //     //     std::cout<<"===="<<std::endl;
    //     // }

    //     // for(int i=0;i<traces_in_span[0].valleys.size();i++)
    //     // {
    //     //     std::cout<<traces_in_span[0].valley_trace_index[i]<<" "<< traces_in_span[0].valleys[i]<<" ";
    //     // }
    //     // std::cout<<std::endl;

    //     std::vector<std::vector<VectorX<T>>> points;
    //     std::vector<std::vector<VectorX<T>>> value;

    //     int num=0;
    //     int pseudo_num=0;
    //     for(int i=0;i<traces_in_span.size();i++)
    //     {
    //         if(!traces_in_span[i].ridges.empty())
    //             num++;
    //         if(!traces_in_span[i].pseudo_ridges.empty())
    //             pseudo_num++;
    //     }

        
    //     points.resize(num);
    //     value.resize(num);


    //     num=0;
    //     for(auto i=0;i<traces_in_span.size();++i)
    //     {
    //         if(traces_in_span[i].ridges.empty())
    //         {
    //             continue;
    //         }
  
    //         points[num].resize(traces_in_span[i].ridges.size());
    //         value[num].resize(traces_in_span[i].ridges.size());


    //         for(auto j=0;j<traces_in_span[i].ridges.size();++j)
    //         {

    //             points[num][j]=core_min+traces_in_span[i].traces[traces_in_span[i].ridge_trace_index[j]][traces_in_span[i].ridges[j]].cwiseProduct(range);
    //             value[num][j].resize(1);
    //             value[num][j].data()[0]=traces_in_span[i].traces_values[traces_in_span[i].ridge_trace_index[j]][traces_in_span[i].ridges[j]];
    //         }

    //         num++;
    //     }

    //     if(num>0)
    //         transfer_data::transfer(points,value,ridge);

    //     points.resize(pseudo_num);
    //     value.resize(pseudo_num);
    //     pseudo_num=0;
    //     for(auto i=0;i<traces_in_span.size();++i)
    //     {
    //         if(traces_in_span[i].pseudo_ridges.empty())
    //         {
    //             continue;
    //         }

    //         points[pseudo_num].resize(traces_in_span[i].pseudo_ridges.size());
    //         value[pseudo_num].resize(traces_in_span[i].pseudo_ridges.size());

    //         for(auto j=0;j<traces_in_span[i].pseudo_ridges.size();++j)
    //         {
    //             points[pseudo_num][j]=core_min+traces_in_span[i].traces[traces_in_span[i].pseudo_ridge_trace_index[j]][traces_in_span[i].pseudo_ridges[j]].cwiseProduct(range);
    //             value[pseudo_num][j].resize(1);
    //             value[pseudo_num][j].data()[0]=traces_in_span[i].traces_values[traces_in_span[i].pseudo_ridge_trace_index[j]][traces_in_span[i].pseudo_ridges[j]];
    //         }

    //         pseudo_num++;
    //     }

    //     point_size=0;
    //     for(auto i = 0;i<traces_in_span.size();++i)
    //     {
    //         point_size += traces_in_span[i].pseudo_ridges.size();
    //     }

    //     if(pseudo_num>0)
    //         transfer_data::transfer(points,value,pseudo_ridge);


    //     if(split_ridge_valley)
    //     {
    //         //valley
    //         point_size=0;
    //         for(auto i = 0;i<traces_in_span.size();++i)
    //         {
    //             point_size += traces_in_span[i].valleys.size();
    //         }

    //         num=0;
    //         for(int i=0;i<traces_in_span.size();i++)
    //         {
    //             if(!traces_in_span[i].valleys.empty())
    //                 num++;
    //         }
            
    //         points.resize(num);
    //         value.resize(num);

    //         num=0;
    //         for(auto i=0;i<traces_in_span.size();++i)
    //         {
    //             if(traces_in_span[i].valleys.empty())
    //             {
    //                 continue;
    //             }

    //             points[num].resize(traces_in_span[i].valleys.size());
    //             value[num].resize(traces_in_span[i].valleys.size());

    //             // auto valley_index = traces_in_span[i].valleys.data();
    //             // auto valley_trace_index=traces_in_span[i].valley_trace_index.data();
    //             // auto point_address = traces_in_span[i].traces.data();
    //             // auto value_address = traces_in_span[i].traces_values.data();

    //             for(auto j=0;j<traces_in_span[i].valleys.size();++j)
    //             {
    //                 points[num][j]= core_min+traces_in_span[i].traces[traces_in_span[i].valley_trace_index[j]][traces_in_span[i].valleys[j]].cwiseProduct(range);
    //                 value[num][j].resize(1);
    //                 value[num][j].data()[0]=traces_in_span[i].traces_values[traces_in_span[i].valley_trace_index[j]][traces_in_span[i].valleys[j]];
    //             }

    //             num++;
    //         }
    //         if(num>0)
    //             transfer_data::transfer(points,value,valley);

    //         //pseudo valley
    //         point_size=0;
    //         for(auto i = 0;i<traces_in_span.size();++i)
    //         {
    //             point_size += traces_in_span[i].pseudo_valleys.size();
    //         }

    //         num=0;
    //         for(int i=0;i<traces_in_span.size();i++)
    //         {
    //             if(!traces_in_span[i].pseudo_valleys.empty())
    //                 num++;
    //         }
            
    //         points.resize(num);
    //         value.resize(num);

    //         num=0;
    //         for(auto i=0;i<traces_in_span.size();++i)
    //         {
    //             if(traces_in_span[i].pseudo_valleys.empty())
    //             {
    //                 continue;
    //             }

    //             points[num].resize(traces_in_span[i].pseudo_valleys.size());
    //             value[num].resize(traces_in_span[i].pseudo_valleys.size());

    //             for(auto j=0;j<traces_in_span[i].pseudo_valleys.size();++j)
    //             {
    //                 points[num][j]= core_min+traces_in_span[i].traces[traces_in_span[i].pseudo_valley_trace_index[j]][traces_in_span[i].pseudo_valleys[j]].cwiseProduct(range);
    //                 value[num][j].resize(1);
    //                 value[num][j].data()[0]=traces_in_span[i].traces_values[traces_in_span[i].pseudo_valley_trace_index[j]][traces_in_span[i].pseudo_valleys[j]];
    //             }

    //             num++;
    //         }
    //         if(num>0)
    //             transfer_data::transfer(points,value,pseudo_valley);

    //     }

  
    // }



template<typename T>
void get_js_error(const Block<T>* b,const Block<T>* b2,
std::vector<TraceInSpan<T>>& traces_in_span)
{
    tbb::affinity_partitioner ap;

    std::vector<std::vector<T>> thread_error(traces_in_span.size());
    std::vector<T> error;

    // tbb::parallel_for(tbb::blocked_range<size_t>(0, traces_in_span.size()),
    // [&](const tbb::blocked_range<size_t>& range)
    // {

        for(auto i=0;i!=traces_in_span.size();++i)
        {
            for(auto j=0;j<traces_in_span[i].traces.size();++j)
            {
                for(auto k=0;k<traces_in_span[i].traces[j].size();++k)
                {
                    T current_error;
                    compute_h(traces_in_span[i].traces[j][k],b,b2,current_error);
                    thread_error[i].emplace_back(std::abs(current_error));
                }
            }
        }
    //  },ap
    // );     



    for(auto i=0;i<traces_in_span.size();++i)
    {
        error.insert(error.end(),thread_error[i].begin(),thread_error[i].end());
    }

    T min_error = *std::min_element(error.begin(), error.end());
    T max_error = *std::max_element(error.begin(), error.end());
    T sum_error = std::accumulate(error.begin(), error.end(), 0.0);
    T average_error = sum_error / error.size();

    printf("min_error: %e, max_error: %e, average_error: %e\n", min_error, max_error, average_error);


}




}

