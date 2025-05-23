#pragma once
#include"trace_in_span.h"

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>
#include <mfa/mfa.hpp>

#include "opts.h"

#include "utility_function.h"

#include "block.hpp"


namespace connect_rv_graph
{
    template<typename T>
    struct EdgeCandidate 
    {
        EdgeCandidate() {}
        EdgeCandidate(int trace_node, T distance, Eigen::VectorX<T> direction)
            : trace_node(trace_node), distance(distance), direction(direction) {}

        int trace_node;             // index of the candidate point 
        T distance;            // distance from the constant point
        Eigen::VectorX<T> direction;  // normalized direction vector from the constant point to the candidate
    };

    template<typename T>
    struct node_info {
        std::vector<EdgeCandidate<T>> edge_candidates;
    };



    template<typename T>
    size_t assign_span_index(const std::vector<std::vector<T>>& span_division, VectorX<T> p, const VectorXi& degree, VectorXi& number_in_every_domain)
    {
        const auto& x_divisions = span_division[0];
        const auto& y_divisions = span_division[1];

        auto x_it = std::lower_bound(x_divisions.begin()+degree[0], x_divisions.end()-degree[0], p[0]);
        auto x_index = std::distance(x_divisions.begin(), x_it) - 1 - degree[0];

        auto y_it = std::lower_bound(y_divisions.begin()+degree[1], y_divisions.end()-degree[1], p[1]);
        auto y_index = std::distance(y_divisions.begin(), y_it) - 1 - degree[1];

        if (x_index < 0 || x_index >= static_cast<int>(x_divisions.size())-2*degree[0]-1 ||
                y_index < 0 || y_index >= static_cast<int>(y_divisions.size()) - 2*degree[0]-1 ) {
                throw std::out_of_range("Point is outside the grid domain");
            }
        
        VectorXi span_index(2);
        span_index << x_index, y_index;


        return utility::obtain_index_from_domain_index(span_index, number_in_every_domain);
        
    }


    //use vector<vector> to store points in every span correspondingly
    template<typename T>
    void reorganize_points(const std::vector<std::vector<T>>& span_division, std::vector<VectorX<T>>& points, std::vector<std::vector<VectorX<T>>>& points_in_span, VectorXi& number_in_every_domain, size_t valid_span_num, size_t real_span_num, const VectorXi& degree,std::unordered_map<size_t, size_t>& real_to_valid_span_index,
    std::vector<std::vector<T>>& critical_point_value,
    const VectorX<T>&         domain_min,
    const VectorX<T>&         domain_range,
    std::vector<size_t>& valid_span_index)
    {
        points_in_span.reserve(real_span_num);
        points_in_span.resize(valid_span_num);
        critical_point_value.reserve(real_span_num);
        critical_point_value.resize(valid_span_num);

        size_t index;
        size_t current_span_num = valid_span_num;
        VectorX<T> point_2d;
        for(auto i=0;i<points.size();++i)
        {
            point_2d = points[i].segment(0,2);
            VectorX<T> point_2d_domain = (point_2d-domain_min).cwiseQuotient(domain_range);
            index = assign_span_index(span_division,point_2d_domain,degree, number_in_every_domain);
            auto it = real_to_valid_span_index.find(index);
            if (it == real_to_valid_span_index.end()) {
                real_to_valid_span_index[index] = current_span_num;

                current_span_num+=1;
                points_in_span.resize(current_span_num);
                points_in_span[current_span_num-1].emplace_back(point_2d);
                critical_point_value[current_span_num-1].emplace_back(points[i].data()[2]);
                valid_span_index.emplace_back(index);
            }
            else
            {
                points_in_span[it->second].emplace_back(point_2d);
                critical_point_value[it->second].emplace_back(points[i].data()[2]);
            }
        }       

    }


    void find_around_span(size_t span_index,VectorXi& number_in_every_domain, const VectorXi& span_num, std::vector<size_t>& around_span_index,std::unordered_map<size_t, size_t>& real_to_valid_span_index)
    {
        VectorXi span_index_2d(2);
        utility::obtainDomainIndex(span_index,span_index_2d,number_in_every_domain);
        VectorXi span_index_2d_around(2);

       
        size_t new_span_index;




        for (int i = -1; i < 2; i++)
        {
            if(span_index_2d[1]+i<0 || span_index_2d[1]+i>=span_num[1])
            {
                continue;
            }
            for(int j=-1;j < 2;j++)
            {
                if(span_index_2d[0]+j<0 || span_index_2d[0]+j>=span_num[0])
                {
                    continue;
                }
                span_index_2d_around << span_index_2d[0]+j,span_index_2d[1]+i;
                new_span_index = utility::obtain_index_from_domain_index(span_index_2d_around,number_in_every_domain);
                auto it = real_to_valid_span_index.find(new_span_index);
                if (it != real_to_valid_span_index.end()) {
                    around_span_index.emplace_back(it->second);               
                }



            }
        }
        
        
        
    }

    template<typename T>
    void connect_critical_points(std::vector<std::vector<VectorX<T>>>& traces,
    std::vector<bool>&is_loop, std::vector<VectorX<T>>& critical_points, int critical_point_span_index, std::vector<std::array<int,8>>& connect_info, T threshold_square,std::vector<T>& min_distance)
    {

        for(int i=0;i<traces.size();i++)
        {
            if(!is_loop[i])
            {
                for(int j=0;j<critical_points.size();j++)
                {
                    T distance[2];
                    distance[0]=(traces[i][0]-critical_points[j]).squaredNorm(); 
                    distance[1]=(traces[i].back()-critical_points[j]).squaredNorm();

                    if(traces[i].size()<4)
                    {
                        int min_index= 0;
                        if(distance[1]<distance[0])
                        {
                            min_index = 1;
                        }

                        // VectorXd example(2);
                        // example << 0.48727,0.0209282;

                        // if((traces[i][0]-example).norm()<1e-4)
                        // {
                        //     std::cout<<"find the point"<<std::endl;
                        //     std::cout<<traces[i].size()<<" "<<critical_points.size()<<std::endl;
                        //     std::cout<<min_index<<std::endl;
                        //     std::cout<<traces[i][0].transpose()<<std::endl;
                        //     std::cout<<distance[min_index]<<" "<<threshold_square <<std::endl;
                        //     std::cout<<min_distance[i]<<" "<<min_distance[i+traces.size()]<<std::endl;
                        //     std::cout<<critical_points[j].transpose()<<std::endl;

                        // }

                        if(distance[min_index]<threshold_square && distance[min_index]<min_distance[i+min_index*traces.size()])
                        {
                            min_distance[i+min_index*traces.size()] = distance[min_index];
                            connect_info[i][min_index*4]=0;
                            connect_info[i][min_index*4+1]=critical_point_span_index;
                            connect_info[i][min_index*4+2]=j;
                        }
                    }
                    else
                    {
                        if(distance[0]<threshold_square && distance[0]<min_distance[i])
                        {
                            min_distance[i] = distance[0];
                            connect_info[i][0]=0;
                            connect_info[i][1]=critical_point_span_index;
                            connect_info[i][2]=j;
                        }
                        if(distance[1]<threshold_square && distance[1]<min_distance[i+traces.size()])
                        {
                            min_distance[i+traces.size()] = distance[1];
                            connect_info[i][4]=0;
                            connect_info[i][5]=critical_point_span_index;
                            connect_info[i][6]=j;
                        }
                    }
                }
            }
        }
    }




    template<typename T>
    void connect_closest_traces(std::vector<size_t>& around_span_index, VectorX<T>& end_point,int* connect_info, std::vector<TraceInSpan<T>>& all_spans, size_t current_span_index, int current_trace_index, T threshold_square, bool test)
    {
        T distance = std::numeric_limits<T>::max();
        T current_distance = distance;

        for(auto i=around_span_index.begin();i<around_span_index.end();++i)
        {
            if(*i>=all_spans.size() || *i==current_span_index)
            {
                continue;
            }
            for(auto j=0;j<all_spans[*i].traces.size();++j)
            {
                if(!all_spans[*i].is_loop[j])
                {
                    if(all_spans[*i].connect_info[j][0]==-1)
                    {
                        current_distance = (end_point-all_spans[*i].traces[j][0]).squaredNorm();
                        if(current_distance<distance)
                        {
                            connect_info[1]=*i;
                            connect_info[2]=j;
                            connect_info[3]=0;
                            distance = current_distance;
                        }
                    }
                    if(all_spans[*i].connect_info[j][4]==-1)
                    {
                        current_distance = (end_point-all_spans[*i].traces[j].back()).squaredNorm();
                        if(current_distance<distance)
                        {
                            connect_info[1]=*i;
                            connect_info[2]=j;
                            connect_info[3]=1;
                            distance = current_distance;
                        }
                    }
                }
            }
        }    

        for(auto j=0;j<all_spans[current_span_index].traces.size();++j)
        {
            if(j!=current_trace_index && (!all_spans[current_span_index].is_loop[j]))
            {
                if(all_spans[current_span_index].connect_info[j][0]==-1)
                {
                    current_distance = (end_point-all_spans[current_span_index].traces[j][0]).squaredNorm();
                    if(current_distance<distance)
                    {
                        connect_info[1]=current_span_index;
                        connect_info[2]=j;
                        connect_info[3]=0;
                        distance = current_distance;
                    }
                }
                if(all_spans[current_span_index].connect_info[j][4]==-1)
                {
                    current_distance = (end_point-all_spans[current_span_index].traces[j].back()).squaredNorm();
                    if(current_distance<distance)
                    {
                        connect_info[1]=current_span_index;
                        connect_info[2]=j;
                        connect_info[3]=1;
                        distance = current_distance;
                    }
                }
            }
        }

        if(test)
        {
            std::cout<<"distance "<<distance<<std::endl;
        }

        if(distance>threshold_square)
        {
            connect_info[1]=-1;
            connect_info[2]=-1;
            connect_info[3]=-1;
        }       

    }



    template<typename T>
    void connect_closest_traces(std::vector<size_t>& around_span_index, TraceInSpan<T>& current_span, std::vector<TraceInSpan<T>>& all_spans, size_t current_span_index, T threshold_square)
    {
        for(auto i=0;i<current_span.traces.size();++i)
        {
            if(!current_span.is_loop[i])
            {
                bool test_start = false;
                // if(std::abs(current_span.traces[i][0][0]-646.296)<1e-2 && std::abs(current_span.traces[i][0][1]-295.531)<1e-2)
                // {
                //     std::cout<<"find the point"<<std::endl;
                //     test_start = true;
                // }

                // if(std::abs(current_span.traces[i].back()[0]-646.296)<1e-2 && std::abs(current_span.traces[i].back()[1]-295.531)<1e-2)
                // {
                //     std::cout<<"find the point 2"<<std::endl;
                //     test_start = true;
                //     std::cout<<current_span.connect_info[i][0]<<" "<<current_span.connect_info[i][1]<<" "<<current_span.connect_info[i][2]<<" "<<current_span.connect_info[i][3]<<std::endl;
                //     std::cout<<current_span.connect_info[i][4]<<" "<<current_span.connect_info[i][5]<<" "<<current_span.connect_info[i][6]<<" "<<current_span.connect_info[i][7]<<std::endl;
                // }

                if(current_span.connect_info[i][0]==-1)
                {
                    connect_closest_traces(around_span_index,current_span.traces[i][0],current_span.connect_info[i].data(),all_spans,current_span_index,i,threshold_square,test_start);
                }

                if(current_span.connect_info[i][4]==-1)
                {
                    connect_closest_traces(around_span_index,current_span.traces[i].back(),current_span.connect_info[i].data()+4,all_spans,current_span_index,i,threshold_square,test_start);
                }
            }
        }

    }


    template<typename T>
    bool angle_near_zero(VectorX<T>& point_need_connect, VectorX<T>& point_already_connect, VectorX<T>& point_to_connect)
    {
        VectorX<T> v1 = point_need_connect-point_already_connect;
        VectorX<T> v2 = point_need_connect-point_to_connect;
        T cos_theta = v1.dot(v2)/(v1.norm()*v2.norm());
        if(cos_theta>0.996)
        {
           return true;
        }
        return false;
    }


 

    template<typename T>
    bool update_connect_info(T point_distance, int min_index, int i, int j, int other_span_index, std::array<int,8>& connect_info, std::vector<T>& max_distance, T threshold_square,int trace_num)
    {
        if(point_distance<threshold_square)
        {
            int start_index = 4*(min_index>>1);
            if(point_distance<max_distance[i+(start_index>>2)*trace_num])
            {
                max_distance[i+(start_index>>2)*trace_num] = point_distance;
                connect_info[start_index]=1;
                connect_info[1+start_index]=other_span_index;
                connect_info[2+start_index]=j;
                connect_info[3+start_index]=min_index%2;
            }
            return true;
        }   
        return false;     
    }


    template<typename T>
    T compute_angle(VectorX<T>& p0, VectorX<T>& p1, VectorX<T>& p2)
    //angle between (p1-p0) and (p2-p0)
    {
        T cos_theta = (p1-p0).normalized().dot((p2-p0).normalized());
        return cos_theta;
    }



    template<typename T>
    bool update_connect_info_for_single_point(T point_distance, int min_index, int i, int j, int other_span_index, std::array<int,8>& connect_info, std::vector<T>& max_distance, T threshold_square,int trace_num,std::vector<std::vector<VectorX<T>>>& critical_points, VectorX<T>& current_point, VectorX<T>& other_point, T theta_threshold)
    {
        if(point_distance<threshold_square)
        {
            if(connect_info[0]!=0 && connect_info[4]!=0)
            {
                if(point_distance<max_distance[i])
                {
                    max_distance[i+trace_num] = max_distance[i];
                    max_distance[i] = point_distance;
                    std::copy(connect_info.begin(), connect_info.begin() + 4, connect_info.begin() + 4);
                    connect_info[0]=1;
                    connect_info[1]=other_span_index;
                    connect_info[2]=j;
                    connect_info[3]=min_index%2;
                }
                else if(point_distance<max_distance[i+trace_num])
                {
                    max_distance[i+trace_num] = point_distance;
                    connect_info[4]=1;
                    connect_info[5]=other_span_index;
                    connect_info[6]=j;
                    connect_info[7]=min_index%2;
                }
                return true;
            }
            else
            {
                if(connect_info[4]!=0)
                {
                    if(point_distance<max_distance[i+trace_num])
                    {
                        // the angle with the edge connecting the critical points cannot be too small
                        VectorX<T>& point = critical_points[connect_info[1]][connect_info[2]];
                        T angle = compute_angle(current_point,point,other_point);
                        if(angle<theta_threshold)
                        {
                            max_distance[i+trace_num] = point_distance;
                            connect_info[4]=1;
                            connect_info[5]=other_span_index;
                            connect_info[6]=j;
                            connect_info[7]=min_index%2;
                        }

                       
                    }
                }
                else if(connect_info[0]!=0)
                {
                    if(point_distance<max_distance[i])
                    {
                        VectorX<T>& point = critical_points[connect_info[5]][connect_info[6]];
                        T angle = compute_angle(current_point,point,other_point);
                        if(angle<theta_threshold)
                        {
                            max_distance[i] = point_distance;
                            connect_info[0]=1;
                            connect_info[1]=other_span_index;
                            connect_info[2]=j;
                            connect_info[3]=min_index%2;
                        }
                    }
                }
                return true;
            }
        }   
        return false;     
    }



    template<typename T>
    void record_connect_info(int case_,std::vector<VectorX<T>>& traces,std::vector<VectorX<T>>& traces_in_other_span, int i, int j, std::vector<T>& max_distance, std::array<int,8>& connect_info, std::array<int,8>& other_connect_info, int other_span_index, T threshold_square,int trace_num,std::vector<std::vector<VectorX<T>>>& critical_points, T theta_threshold)
    {


        switch (case_)
        {
        case 4: //if the trace only contain one point, need to consider it seperately 
        {
            T point_distance[2];
            point_distance[0] = (traces[0]-traces_in_other_span[0]).squaredNorm();
            point_distance[1] = (traces[0]-traces_in_other_span.back()).squaredNorm();
            int min_index = 0;
            if(point_distance[1]<point_distance[0])
            {
                min_index = 1;
            }

           

            update_connect_info_for_single_point(point_distance[min_index],min_index,i,j,other_span_index,connect_info,max_distance,threshold_square,trace_num,critical_points,traces[0],traces_in_other_span[min_index*(traces_in_other_span.size()-1)],theta_threshold);

            if(traces_in_other_span.size()>=4)
            {
                int corresponding_index = 1-min_index;
                update_connect_info_for_single_point(point_distance[corresponding_index],corresponding_index,i,j,other_span_index,other_connect_info,max_distance,threshold_square,trace_num,critical_points,traces[0],traces_in_other_span[corresponding_index*(traces_in_other_span.size()-1)],theta_threshold);
            }

        }
            break;

        case 3: // both traces don't connect to critical points
        { 
            if(traces.size()<4 && traces_in_other_span.size()<4)
            {                              
        
                T point_distance[4];
                point_distance[0] = (traces[0]-traces_in_other_span[0]).squaredNorm();
                point_distance[1] = (traces[0]-traces_in_other_span.back()).squaredNorm();
                point_distance[2] = (traces.back()-traces_in_other_span[0]).squaredNorm();
                point_distance[3] = (traces.back()-traces_in_other_span.back()).squaredNorm();

                int min_index = std::distance(std::begin(point_distance),
                    std::min_element(std::begin(point_distance), std::end(point_distance)));
                
                bool test = update_connect_info(point_distance[min_index],min_index,i,j,other_span_index,connect_info,max_distance,threshold_square,trace_num);  

            //      VectorXd example2(2);
            // example2 << 227.441, 152.891;
            // VectorXd example(2);
            // example << 227.287, 151.462;


            // if((traces[0]-example).norm()<1e-2)
            // {
            //     std::cout<<"find the point"<<std::endl;
            //     std::cout<<traces.size()<<" "<<traces_in_other_span.size()<<std::endl;
            //     std::cout<<traces_in_other_span[0].transpose()<<" "<< traces_in_other_span.back().transpose()<<std::endl;
            //     std::cout<<min_index<<std::endl;
            //     std::cout<<traces[0].transpose()<<std::endl;
            //     std::cout<<point_distance[min_index]<<" "<<threshold_square <<std::endl;

            // }

            }
            else
            {
               T point_distance[4];
                point_distance[0] = (traces[0]-traces_in_other_span[0]).squaredNorm();
                point_distance[1] = (traces[0]-traces_in_other_span.back()).squaredNorm();
                point_distance[2] = (traces.back()-traces_in_other_span[0]).squaredNorm();
                point_distance[3] = (traces.back()-traces_in_other_span.back()).squaredNorm();

                int min_index = std::distance(std::begin(point_distance),
                    std::min_element(std::begin(point_distance), std::end(point_distance)));


                
            // VectorXd example2(2);
            // example2 << 227.441, 152.891;
            // VectorXd example(2);
            // example << 227.287, 151.462;


            // if((traces[0]-example).norm()<1e-2)
            // {
            //     std::cout<<"find the point"<<std::endl;
            //     std::cout<<traces.size()<<" "<<traces_in_other_span.size()<<std::endl;
            //     std::cout<<traces_in_other_span[0].transpose()<<" "<< traces_in_other_span.back().transpose()<<std::endl;
            //     std::cout<<min_index<<std::endl;
            //     std::cout<<traces[0].transpose()<<std::endl;
            //     std::cout<<point_distance[min_index]<<" "<<threshold_square <<std::endl;
    
            // }
                
                bool test_1=update_connect_info(point_distance[min_index],min_index,i,j,other_span_index,connect_info,max_distance,threshold_square,trace_num);
               
                int corresponding_index = 3-min_index;
                bool test_2 = update_connect_info(point_distance[corresponding_index],corresponding_index,i,j,other_span_index,connect_info,max_distance,threshold_square,trace_num);




            }

        }
            break;
        
        case 1: // this trace has one critical point
        {
            int avaliable_index = 0; // 0: start, 1: end
            if(connect_info[0]==0)
            {
                avaliable_index = traces.size()-1;
            }            
            T point_distance[2];
            point_distance[0] = (traces[avaliable_index]-traces_in_other_span[0]).squaredNorm();
            point_distance[1] = (traces[avaliable_index]-traces_in_other_span.back()).squaredNorm();
            int min_index = 0;
            if(point_distance[1]<point_distance[0])
            {
                min_index = 1;
            }
            int compatible_index = min_index+2*(avaliable_index/(traces.size()-1));
            bool test = update_connect_info(point_distance[min_index],compatible_index,i,j,other_span_index,connect_info,max_distance,threshold_square,trace_num);

        }
            break;
        case 2: // other trace has one critical point
        {
            int avaliable_index = 0; // 0: start, 1: end
            if(other_connect_info[0]==0)
            {
                avaliable_index = traces_in_other_span.size()-1;
            }
            T point_distance[2];
            point_distance[0] = (traces[0]-traces_in_other_span[avaliable_index]).squaredNorm();
            point_distance[1] = (traces.back()-traces_in_other_span[avaliable_index]).squaredNorm();
            int min_index = 0;
            if(point_distance[1]<point_distance[0])
            {
                min_index = 1;
            }
            int compatible_index;

            if(traces_in_other_span.size()==1)
            {
                compatible_index = 2*min_index;
            }
            else
            {
                compatible_index = 2*min_index+avaliable_index/(traces_in_other_span.size()-1);
            }

            bool test = update_connect_info(point_distance[min_index],compatible_index,i,j,other_span_index,connect_info,max_distance,threshold_square,trace_num);

           
        }
            break;
        case 0: //  both of them have one critical point
        {
            int avaliable_index = 0; // 0: start, 1: end
            if(connect_info[0]==0)
            {
                avaliable_index = traces.size()-1;
            }
            int other_avaliable_index = 0; // 0: start, 1: end
            if(other_connect_info[0]==0)
            {
                other_avaliable_index = traces_in_other_span.size()-1;
            }
            T point_distance = (traces[avaliable_index]-traces_in_other_span[other_avaliable_index]).squaredNorm();

            int compatible_index;
            if(traces_in_other_span.size()==1)
            {
                compatible_index = 2*avaliable_index/(traces.size()-1)+1;
            }
            else{
                compatible_index = 2*avaliable_index/(traces.size()-1) + other_avaliable_index/(traces_in_other_span.size()-1);
            }

            bool test = update_connect_info(point_distance,compatible_index,i,j,other_span_index,connect_info,max_distance,threshold_square,trace_num);

            // VectorXd example(2);
            // example << 0.48727,0.0209282;
            // if((traces[0]-example).norm()<1e-4)
            // {
            //     std::cout<<"find the point +++"<<std::endl;
            //     std::cout<<traces.size()<<" "<<traces_in_other_span.size()<<std::endl;
            //     std::cout<<traces[0].transpose()<<std::endl;
            //     std::cout<<compatible_index<<" "<<avaliable_index<<std::endl;
            //     std::cout<<connect_info[0]<<" "<<connect_info[1]<<" "<<connect_info[2]<<" "<<connect_info[3]<<" "<<connect_info[4]<<" "<<connect_info[5]<<" "<<connect_info[6]<<" "<<connect_info[7]<<std::endl;

            // }

        }
            break;
        }
    }





    template<typename T>
    void connect_other_traces(std::vector<std::vector<VectorX<T>>>& traces,std::vector<bool>&is_loop, std::vector<std::vector<VectorX<T>>>& traces_in_other_span, std::vector<bool>&other_is_loop, int other_span_index, std::vector<std::array<int,8>>& connect_info,std::vector<std::array<int,8>>& other_span_connect_info, int current_span_index, T threshold_square, std::vector<T>& max_distance,std::vector<size_t>& around_span_index,std::vector<std::vector<VectorX<T>>>& critical_points, T theta_threshold)
    {

        for(int i=0;i<traces.size();i++)
        {
            if(!is_loop[i])
            {
                int trace_connect_index = 0; // 1: both not connected. 0. one node is connnected
                if(connect_info[i][0]==0 && connect_info[i][4]==0)
                {
                    continue;
                }

                if(connect_info[i][0]!=0 && connect_info[i][4]!=0)
                {
                    trace_connect_index = 1;
                }



                for(int j = 0; j < traces_in_other_span.size(); j++)
                {
                    

                    if(!other_is_loop[j])
                    {

                        if(other_span_connect_info[j][0]==0 && other_span_connect_info[j][4]==0)
                        {
                            continue;
                        }

                        int other_trace_connect_index = 0; // 1: both not connected. 0. one node is connnected
                        if(other_span_connect_info[j][0]!=0 && other_span_connect_info[j][4]!=0)
                        {
                            other_trace_connect_index = 1;
                        }

                        int case_ = 2*trace_connect_index+other_trace_connect_index;
                        if(traces[i].size()==1)
                        {
                            case_ = 4;
                        }       

                        // VectorXd example2(2);
                        // example2 << -0.39243, 0.0341784;
                        // if((traces[i].back()-example2).norm()<1e-4)
                        // {
                        //     std::cout<<"find the point"<<std::endl;
                        //     std::cout<< traces[i].back().transpose()<<std::endl;
                        //     std::cout<<traces_in_other_span[j].back().transpose()<<std::endl;
                        //     std::cout<<traces_in_other_span[j][0].transpose()<<std::endl;
                        //     std::cout<<traces[i].size()<<" "<<traces_in_other_span[j].size()<<std::endl;
                        //     std::cout<<trace_connect_index<<" "<<other_trace_connect_index<<" "<<case_<<std::endl;
                        // }

     
                        //case_: 0: both traces have one avalible node, 1: this trace has one avalible node, 2: other trace has one avalible node, 3: both of them are not connected, 4: the trace only has one node
                        record_connect_info(case_,traces[i],traces_in_other_span[j],i,j,max_distance,connect_info[i],other_span_connect_info[j],other_span_index,threshold_square,traces.size(),critical_points,theta_threshold);
                    
                    }

                }
                
            }
        }

    }


    template<typename T>
    void connect_inside_span(std::vector<std::vector<VectorX<T>>>& traces, std::vector<bool>&is_loop, std::vector<std::array<int,8>>& connect_info, int current_span_index,T threshold_square,std::vector<T>& max_distance,T step_size_square,std::vector<std::vector<VectorX<T>>>& critical_points, T theta_threshold)
    {

        for(int i=0;i<traces.size();i++)
        {
            if(!is_loop[i])
            {
                int trace_connect_index = 0; // 1: both not connected. 0. one node is connnected
                if(connect_info[i][0]==0 && connect_info[i][4]==0)
                {
                    continue;
                }

                if(connect_info[i][0]!=0 && connect_info[i][4]!=0)
                {
                    trace_connect_index = 1;
                }

                for(int j = 0; j < traces.size(); j++)
                {

                    if(j!=i && (!is_loop[j]))
                    {

                        if(connect_info[j][0]==0 && connect_info[j][4]==0)
                        {
                            continue;
                        }

                        int other_trace_connect_index = 0; // 1: both not connected. 0. one node is connnected
                        if(connect_info[j][0]!=0 && connect_info[j][4]!=0)
                        {
                            other_trace_connect_index = 1;
                        }

                        int case_ = 2*trace_connect_index+other_trace_connect_index;
                        if(traces[i].size()==1)
                        {
                            case_ = 4;
                        }
                        //case_: 0: both traces have one avalible node, 1: this trace has one avalible node, 2: other trace has one avalible node, 3: both of them are not connected, 4: the trace only has one node
                        record_connect_info(case_,traces[i],traces[j],i,j,max_distance,connect_info[i],connect_info[j],current_span_index,threshold_square,traces.size(),critical_points,theta_threshold);


                    }

                }
                
            }
        }
    }





    template<typename T>
    void compare_with_closed_critical_point(std::vector<size_t>& around_span_index, VectorX<T>& end_point, int* connect_info,std::vector<std::vector<VectorX<T>>>& critical_points, T threshold_square)
    {
        T distance = std::numeric_limits<T>::max();

        T current_distance = distance;

        for(auto i=around_span_index.begin();i<around_span_index.end();++i)
        {
            if(!critical_points[*i].empty())
            {
                for(auto j=0;j<critical_points[*i].size();++j)
                {
                    current_distance = (end_point-critical_points[*i][j]).squaredNorm();
                    if(current_distance<distance)
                    {
                        connect_info[1]=*i;
                        connect_info[2]=j;
                        distance = current_distance;
                    }
                }
            }
        }

        if(distance<threshold_square)
        {
            connect_info[0]=0;
        }
        else
        {
            connect_info[1]=-1;
            connect_info[2]=-1;
        }

    }

    template<typename T>
    void compare_with_closed_critical_point(std::vector<std::vector<VectorX<T>>>& critical_points, TraceInSpan<T>& span,std::vector<size_t>& around_span_index,T threshold_square)
    {
        for(auto i=0;i<span.traces.size();++i)
        {
            if(!span.is_loop[i])
            {
                if(span.connect_info[i][0]==-1)
                {
                    compare_with_closed_critical_point(around_span_index,span.traces[i][0],span.connect_info[i].data(),critical_points,threshold_square);
                }

                if(span.connect_info[i][4]==-1)
                {
                    compare_with_closed_critical_point(around_span_index,span.traces[i].back(),span.connect_info[i].data()+4,critical_points,threshold_square);
                }
            }
        }
    }




    template<typename T>
    void compare_with_neighbor_span_critical_point(std::vector<std::vector<VectorX<T>>>& critical_points, TraceInSpan<T>& span, size_t current_valid_span_index, std::vector<TraceInSpan<T>>& all_spans, VectorXi& number_in_every_domain, const VectorXi& span_num,
    T threshold_square, std::unordered_map<size_t, size_t>& real_to_valid_span_index,std::vector<size_t>& valid_span_index,
    std::vector<size_t>& around_span_index,std::vector<size_t>& critical_points_prefix, T cosThreshold)
    {
        // around_span_index.reserve(9);

        // find_around_span(valid_span_index[current_valid_span_index],number_in_every_domain,span_num,around_span_index,real_to_valid_span_index);

        // std::cout<<"test 1"<<std::endl;

        std::vector<T> min_distance(2*span.traces.size(),std::numeric_limits<T>::max());

       


        for(auto i=around_span_index.begin();i<around_span_index.end();++i)
        {
            if(!critical_points[*i].empty())
            {
                VectorXi span_index_2d(2);



                utility::obtainDomainIndex(valid_span_index[*i],span_index_2d,number_in_every_domain);
                VectorXd p0(2);

                // std::cout<<"test 1.1"<<std::endl;

                // std::cout<<span_index_2d[0]+mfa_data->p[0]<<" "<<span_index_2d[1]<<" "<<mfa_data->p[1]<<std::endl;

                // p0[0]=mfa_data->tmesh.all_knots[0][span_index_2d[0]+mfa_data->p[0]];
                // p0[1]=mfa_data->tmesh.all_knots[1][span_index_2d[1]+mfa_data->p[1]];
                // // std::cout<<"test 1.2"<<std::endl;
                // VectorXd p1(2);
                // p1[0]=mfa_data->tmesh.all_knots[0][span_index_2d[0]+1+mfa_data->p[0]];
                // p1[1]=mfa_data->tmesh.all_knots[1][span_index_2d[1]+1+mfa_data->p[1]];

                // // std::cout<<"test 1.3"<<std::endl;
                // VectorXd p0_=p0.cwiseProduct(domain_range)+domain_min;
                // VectorXd p1_=p1.cwiseProduct(domain_range)+domain_min;

                // std::cout<<"test 2"<<std::endl;

                connect_critical_points(span.traces,span.is_loop,critical_points[*i],*i,span.connect_info, threshold_square,min_distance);

                // std::cout<<"test 3"<<std::endl;
            }        
            
        }

        std::unordered_map<size_t, node_info<T>> critical_point_map;

        for(int i=0;i<span.traces.size();i++)
        {
            if(!span.is_loop[i])
            {
                if(span.connect_info[i][0]==0)
                {
                    

                    size_t critical_point_index = critical_points_prefix[span.connect_info[i][1]]+span.connect_info[i][2];
                    auto key=critical_point_map.find(critical_point_index);
                    if(key==critical_point_map.end())
                    {

                        node_info<T> temp_node_info;
                        EdgeCandidate<T> temp_edge_candidate;
                        temp_edge_candidate.trace_node = i;
                        temp_edge_candidate.distance = (span.traces[i][0]-critical_points[span.connect_info[i][1]][span.connect_info[i][2]]).squaredNorm();
                        temp_edge_candidate.direction = (span.traces[i][0]-critical_points[span.connect_info[i][1]][span.connect_info[i][2]]).normalized();
                        temp_node_info.edge_candidates.emplace_back(temp_edge_candidate);
                        critical_point_map[critical_point_index] = temp_node_info;
                    }
                    else
                    {
                        key->second.edge_candidates.emplace_back(EdgeCandidate<T>{i,(span.traces[i][0]-critical_points[span.connect_info[i][1]][span.connect_info[i][2]]).squaredNorm(),(span.traces[i][0]-critical_points[span.connect_info[i][1]][span.connect_info[i][2]]).normalized()});
                    }
                }
                if(span.connect_info[i][4]==0)
                {
                    size_t critical_point_index = critical_points_prefix[span.connect_info[i][5]]+span.connect_info[i][6];
                    auto key=critical_point_map.find(critical_point_index);
                    if(key==critical_point_map.end())
                    {
                        node_info<T> temp_node_info;
                        EdgeCandidate<T> temp_edge_candidate;
                        temp_edge_candidate.trace_node = i+span.traces.size();
                        temp_edge_candidate.distance = (span.traces[i].back()-critical_points[span.connect_info[i][5]][span.connect_info[i][6]]).squaredNorm();
                        temp_edge_candidate.direction = (span.traces[i].back()-critical_points[span.connect_info[i][5]][span.connect_info[i][6]]).normalized();
                        temp_node_info.edge_candidates.emplace_back(temp_edge_candidate);
                        critical_point_map[critical_point_index] = temp_node_info;
                    }
                    else
                    {
                        key->second.edge_candidates.emplace_back(EdgeCandidate<T>{i+(int)(span.traces.size()),(span.traces[i].back()-critical_points[span.connect_info[i][5]][span.connect_info[i][6]]).squaredNorm(),(span.traces[i].back()-critical_points[span.connect_info[i][5]][span.connect_info[i][6]]).normalized()});
                    }
                }
            }
        }

        for(auto& entry : critical_point_map)
        {
            size_t key = entry.first;
            if(entry.second.edge_candidates.size()>1)
            {
                std::vector<EdgeCandidate<T>>& candidate = entry.second.edge_candidates;
                std::vector<bool> overlap(candidate.size(),false);

                std::sort(candidate.begin(), candidate.end(), [](const EdgeCandidate<T> &a, const EdgeCandidate<T> &b) {
                    return a.distance < b.distance;
                });

                for (auto m=1;m<candidate.size();++m) {
                    // Compare candidate's direction with all accepted edges.
                    for (auto k=0;k< m;++k) {
                        if(overlap[k])
                        {
                            continue;
                        }
                        T dot = candidate[m].direction.dot(candidate[k].direction);
                        // If the dot product is higher than cosThreshold, the angle is below the threshold.
                        if (dot > cosThreshold) {
                            overlap[m] = true;
                            break;
                        }
                    }
                }

                for(auto m=1;m<candidate.size();++m)
                {
                    if(overlap[m])
                    {
                        span.connect_info[candidate[m].trace_node % span.traces.size()][4*(candidate[m].trace_node / span.traces.size())] =-1;                      
                    }
                }
            }
        }




    }

    //record 1 to connect_info [0] and [4]
    void complete_connect_info(std::vector<std::array<int,8>>& connect_info)
    {
        for(int i=0;i<connect_info.size();i++)
        {
            if(connect_info[i][0]==-1 && connect_info[i][1]!=-1)
            {
                connect_info[i][0]=1;
            }
            if(connect_info[i][4]==-1 && connect_info[i][5]!=-1)
            {
                connect_info[i][4]=1;
            }
        }
    }



    template<typename T>
    bool find_duplicate(TraceInSpan<T>& span,std::vector<size_t>& around_span_index,T step_size_square,std::vector<TraceInSpan<T>>& all_spans, int i)
    {
        for(auto j=around_span_index.begin();j<around_span_index.end();j++)
        {
            if(*j>=all_spans.size())
            {
                continue;
            }
            for(int k=0;k<all_spans[*j].traces.size();++k)
            {
                
                if(all_spans[*j].is_loop[k])
                {
                    if(utility::check_if_a_trace_is_in_a_loop(all_spans[*j].traces[k],span.traces[i],step_size_square))
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    template<typename T>
    T distanceSquaredToSegment(VectorX<T>& point, VectorX<T>& start, VectorX<T>& end)
    {
        VectorX<T> v = end - start;
        VectorX<T> w = point - start;

        T c1 = w.dot(v);
        if (c1 <= 0) {
            return w.squaredNorm();
        }

        T c2 = v.dot(v);
        if (c2 <= c1) {
            return (point - end).squaredNorm();
        }

        T b = c1 / c2;
        VectorX<T> pb = start + b * v;
        return (point - pb).squaredNorm();
    }


    template<typename T>
    bool check_if_a_point_is_in_a_trace(std::vector<VectorX<T>>& trace, VectorX<T>& point, T threshold_square)
    {
        T distance_threshold = 0.05*0.05*threshold_square;

        int closest_distance = -1;
        for(int i=1;i<trace.size();i++)
        {
            T distance = distanceSquaredToSegment(point,trace[i-1],trace[i]);

            // std::cout<<distance<<" "<<distance_threshold<<trace.size()<<std::endl;
            // std::cout<<trace[i-1].transpose()<<" "<<trace[i].transpose()<<" "<<point.transpose()<<std::endl;
            if(distance<distance_threshold)
            {
                return true;
            }
        }
        return false;
    }




    template<typename T>
    bool find_node_on_edge(TraceInSpan<T>& span,std::vector<size_t>& around_span_index,T step_size_square,std::vector<TraceInSpan<T>>& all_spans, int i)
    {
        for(auto j=around_span_index.begin();j<around_span_index.end();j++)
        {
            if(*j>=all_spans.size())
            {
                continue;
            }
            for(int k=0;k<all_spans[*j].traces.size();++k)
            {
                if((!all_spans[*j].is_loop[k])&&all_spans[*j].traces[k].size()>1)
                {
                    if(utility::check_if_a_trace_is_in_a_loop(all_spans[*j].traces[k],span.traces[i],step_size_square))
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }


    template<typename T>
    //delete the duplicate edges in loops and delete single nodes in an edge
    void delete_duplicate_edges_and_nodes(TraceInSpan<T>& span,std::vector<size_t>& around_span_index,T step_size_square,std::vector<TraceInSpan<T>>& all_spans)
    {
        std::vector<bool> traces_should_remove(span.traces.size(),false);

        

        for(int i=0;i<span.traces.size();++i)
        {
            if(!span.is_loop[i])
            {
                if(find_duplicate(span,around_span_index,step_size_square,all_spans,i))
                {
                    traces_should_remove[i] = true;
                }
                else
                {
                    if(span.traces[i].size()==1)
                    {
                        if(find_node_on_edge(span,around_span_index,step_size_square,all_spans,i))
                        {
                            traces_should_remove[i] = true;
                        }
                        
                    }
                }
            }
        }




        bool all_false = std::all_of(traces_should_remove.begin(), traces_should_remove.end(), [](bool val) {
            return !val;
        });
        
        if(!all_false)
        {
            TraceInSpan<T> temp_traces;
            for (int i = 0; i < span.traces.size(); ++i) {
                if (!traces_should_remove[i]) {
                    temp_traces.traces.emplace_back(span.traces[i]);
                    temp_traces.is_loop.emplace_back(span.is_loop[i]);
                    temp_traces.traces_values.emplace_back(span.traces_values[i]);
                }
            }
            if(!span.point_types.empty())
            {
                for (int i = 0; i < span.traces.size(); ++i) {
                    if (!traces_should_remove[i]) {
                    temp_traces.point_types.emplace_back(span.point_types[i]);
                    }
                }
            }
            if(!temp_traces.traces.empty())
            {
                span = temp_traces;
            }
            else
            {
                span.traces.clear();
                span.is_loop.clear();
                span.traces_values.clear();
                span.point_types.clear();
            }
        }

    }



    template<typename T>
    void connect_trace_two_ends(TraceInSpan<T>& span,T step_size_square,size_t current_span_index)
    {
        
        for(int i=0;i<span.traces.size();++i)
        {
            if(span.connect_info[i][0]==-1 && span.connect_info[i][4]==-1)
            {
                if(span.traces[i].size()>2 && (span.traces[i][0]-span.traces[i].back()).squaredNorm()<step_size_square)
                {
                    span.connect_info[i][0]=1;
                    span.connect_info[i][1]=current_span_index;
                    span.connect_info[i][2]=i;
                    span.connect_info[i][3]=1;
                    span.connect_info[i][4]=1;
                    span.connect_info[i][5]=current_span_index;
                    span.connect_info[i][6]=i;
                    span.connect_info[i][7]=0;
                }
            }
        }
        
    }

    template<typename T>
    void compare_with_neighbor_span_traces(TraceInSpan<T>& span, size_t current_valid_span_index, std::vector<TraceInSpan<T>>& all_spans,
    T threshold_square,
    std::vector<size_t>& around_span_index,T step_size_square,std::vector<std::vector<VectorX<T>>>& critical_points, T theta_threshold)
    {
        std::vector<T> max_distance(2*span.traces.size(),std::numeric_limits<T>::max());


        for(auto i=around_span_index.begin();i<around_span_index.end();++i)
        {
            if(*i==current_valid_span_index || *i>=all_spans.size())
            {
                continue;
            }
            connect_other_traces(span.traces,span.is_loop,all_spans[*i].traces,all_spans[*i].is_loop,*i,span.connect_info,all_spans[*i].connect_info,current_valid_span_index,threshold_square,max_distance,around_span_index,critical_points,theta_threshold);            
        }

        connect_inside_span(span.traces,span.is_loop,span.connect_info,current_valid_span_index,threshold_square,max_distance,step_size_square,critical_points,theta_threshold);
    }





    template<typename T>
    void construct_edge(int* connect_info, std::vector<std::array<size_t,2>>& edge, size_t current_end_index,
    const std::vector<size_t>& critical_points_prefix,const std::vector<TraceInSpan<T>>& all_ridge_valley_in_span,
    std::vector<std::shared_ptr<std::atomic<int>>>& vertex_count)
    {

        if(connect_info[0]==0)
        {                       
            edge.emplace_back(std::array<size_t,2>{current_end_index,critical_points_prefix[connect_info[1]]+connect_info[2]});
            vertex_count[current_end_index]->fetch_add(1);
            vertex_count[critical_points_prefix[connect_info[1]]+connect_info[2]]->fetch_add(1);
        }
        else if(connect_info[0]==1)
        {

            size_t index = all_ridge_valley_in_span[connect_info[1]].prefix_start[connect_info[2]]+
            connect_info[3]*(all_ridge_valley_in_span[connect_info[1]].traces[connect_info[2]].size()-1);


            if(index>current_end_index)
            {

                edge.emplace_back(std::array<size_t,2>{current_end_index,index});   
                vertex_count[current_end_index]->fetch_add(1);
                vertex_count[index]->fetch_add(1);     
            }                               
        }
    }

    template<typename T>
    void construct_edge_test(int* connect_info, std::vector<std::array<size_t,2>>& edge, size_t current_end_index,
    const std::vector<size_t>& critical_points_prefix,const std::vector<TraceInSpan<T>>& all_ridge_valley_in_span,
    std::vector<std::shared_ptr<std::atomic<int>>>& vertex_count, int record_valid_span_index)
    {


        if(connect_info[0]==0)
        {          
            std::cout<<"loop "<<record_valid_span_index<<std::endl;             
            // edge.emplace_back(std::array<size_t,2>{current_end_index,critical_points_prefix[connect_info[1]]+connect_info[2]});
            // vertex_count[current_end_index]->fetch_add(1);
            // vertex_count[critical_points_prefix[connect_info[1]]+connect_info[2]]->fetch_add(1);
        }
        else if(connect_info[0]==1)
        {
            size_t index = all_ridge_valley_in_span[connect_info[1]].prefix_start[connect_info[2]]+
            connect_info[3]*(all_ridge_valley_in_span[connect_info[1]].traces[connect_info[2]].size()-1);

            std::cout<< current_end_index<<" "<<index<<" "<<record_valid_span_index<<std::endl;

            // if(index>current_end_index)
            // {
            //     edge.emplace_back(std::array<size_t,2>{current_end_index,index});   
            //     vertex_count[current_end_index]->fetch_add(1);
            //     vertex_count[index]->fetch_add(1);     
            // }                               
        }

    }

    //supplement for the above function, check index<current_end_index
    template<typename T>
    void construct_edge_sup(int* connect_info, std::vector<std::array<size_t,2>>& edge, size_t current_end_index,
    const std::vector<TraceInSpan<T>>& all_ridge_valley_in_span,
    std::vector<std::shared_ptr<std::atomic<int>>>& vertex_count)
    {
        if(connect_info[0]==1 && vertex_count[current_end_index]->load()<2)
        {
            size_t index = all_ridge_valley_in_span[connect_info[1]].prefix_start[connect_info[2]]+
            connect_info[3]*(all_ridge_valley_in_span[connect_info[1]].traces[connect_info[2]].size()-1);

            std::cout<<"add index 1"<<index<<" current_end_index "<<current_end_index<<std::endl;

            if(index<current_end_index)
            {
                std::cout<<"add index "<<index<<" current_end_index "<<current_end_index<<std::endl;
                edge.emplace_back(std::array<size_t,2>{current_end_index,index});   
                vertex_count[current_end_index]->fetch_add(1);
                vertex_count[index]->fetch_add(1);     
            }
        }
    }


    template<typename T>
    void construct_edge(std::vector<std::array<size_t,2>>& edge, //tbb::concurrent_vector
                TraceInSpan<T>& ridge_valley_in_span,
                std::vector<TraceInSpan<T>>& all_ridge_valley_in_span,
                const std::vector<std::vector<VectorX<T>>>& critical_points_in_span,
                const std::vector<size_t>& critical_points_prefix,
                std::vector<std::shared_ptr<std::atomic<int>>>& vertex_count, int record_valid_span_index)
    {
        for(auto i=0;i<ridge_valley_in_span.traces.size();++i)
        {
            for(int j=1;j<ridge_valley_in_span.traces[i].size();++j)
            {
                edge.emplace_back(std::array<size_t,2>{ridge_valley_in_span.prefix_start[i]+j-1,ridge_valley_in_span.prefix_start[i]+j});   
                vertex_count[ridge_valley_in_span.prefix_start[i]+j-1]->fetch_add(1);
                vertex_count[ridge_valley_in_span.prefix_start[i]+j]->fetch_add(1);             
            }
            if(ridge_valley_in_span.is_loop[i])
            {
                edge.emplace_back(std::array<size_t,2>{ridge_valley_in_span.prefix_start[i]+ridge_valley_in_span.traces[i].size()-1,ridge_valley_in_span.prefix_start[i]});
                vertex_count[ridge_valley_in_span.prefix_start[i]+ridge_valley_in_span.traces[i].size()-1]->fetch_add(1);
                vertex_count[ridge_valley_in_span.prefix_start[i]]->fetch_add(1);
            }
            else //connect with critical points and other traces
            {               


                // VectorXd example(2);
                // example << 0.48727,0.0209282;

                // if((ridge_valley_in_span.traces[i][0]-example).norm()<1e-4)
                // {
                //     std::cout<<"find the point"<<std::endl;
                //     std::cout<<ridge_valley_in_span.traces[i].size()<<std::endl;
                //     std::cout<<ridge_valley_in_span.traces[i][0].transpose()<<std::endl;
                //     std::cout<<ridge_valley_in_span.connect_info[i][0]<<" "<<ridge_valley_in_span.connect_info[i][1]<<" "<<ridge_valley_in_span.connect_info[i][2]<<" "<<ridge_valley_in_span.connect_info[i][3]<<std::endl;

                // }

                //check edge with only 1 vertex
                // left
                auto connect_info = ridge_valley_in_span.connect_info[i].data();
                construct_edge(connect_info,edge,ridge_valley_in_span.prefix_start[i],critical_points_prefix,all_ridge_valley_in_span,vertex_count);
                //right

                if(ridge_valley_in_span.traces[i].size()==1 && connect_info[0]==connect_info[4] && connect_info[1] == connect_info[5] && connect_info[2] == connect_info[6] && connect_info[3] == connect_info[7])
                {
                    continue;
                }
                connect_info = ridge_valley_in_span.connect_info[i].data() + 4;
                construct_edge(connect_info,edge,ridge_valley_in_span.prefix_start[i]+ridge_valley_in_span.traces[i].size()-1,critical_points_prefix,all_ridge_valley_in_span,vertex_count);     

                // if(record_valid_span_index>0)
                // {
                //     construct_edge_test(connect_info,edge,ridge_valley_in_span.prefix_start[i]+ridge_valley_in_span.traces[i].size()-1,critical_points_prefix,all_ridge_valley_in_span,vertex_count,record_valid_span_index);   
                //     connect_info = ridge_valley_in_span.connect_info[i].data(); 
                //     construct_edge_test(connect_info,edge,ridge_valley_in_span.prefix_start[i],critical_points_prefix,all_ridge_valley_in_span,vertex_count,record_valid_span_index);
                // }           
            }           
        }



    }



    template<typename T>
    void construct_edge(std::vector<std::array<size_t,2>>& edge, //tbb::concurrent_vector
                std::vector<TraceInSpan<T>>& ridge_valley_in_span,
                const std::vector<std::vector<VectorX<T>>>& critical_points_in_span,
                const std::vector<size_t>& critical_points_prefix, int vertex_num, int record_valid_span_index, int record_valid_span_index2)
    {
        std::vector<std::shared_ptr<std::atomic<int>>> vertex_count(vertex_num,std::make_shared<std::atomic<int>>(0));

        for(auto i=0;i<ridge_valley_in_span.size();++i)
        {
            int test = 0;
            if(i==record_valid_span_index)
            {
                test = 1;
            }
            if(i==record_valid_span_index2)
            {
                test = 2;
            }
           
            construct_edge(edge,ridge_valley_in_span[i],ridge_valley_in_span,critical_points_in_span,critical_points_prefix,vertex_count,test);
        }     



        // for(auto i=ridge_valley_in_span.begin();i<ridge_valley_in_span.end();++i)
        // {
        //     for(auto j=0;j<i->traces.size();++j)
        //     {
        //         construct_edge_sup(i->connect_info[j].data(),edge,i->prefix_start[j],ridge_valley_in_span,vertex_count);
        //         construct_edge_sup(i->connect_info[j].data()+4,edge,i->prefix_start[j]+i->traces[j].size()-1,ridge_valley_in_span,vertex_count);
        //     }
        // }
    }



    template<typename T>
    void writeEdgeValue(
                const std::string& edge_value_filename,
                const std::vector<TraceInSpan<T>>& ridge_valley_in_span,
                const std::vector<std::vector<VectorX<T>>>& critical_points_in_span,
                const std::vector<std::array<size_t,2>>& edges) 
    {

        
        std::vector<int> vertex_type;
        // Write vertex data
       

        for (auto i=ridge_valley_in_span.begin();i<ridge_valley_in_span.end();++i)
        {
            for(auto j=0;j<i->traces.size();++j)
            {
                for(auto k=0;k<i->traces[j].size();++k)
                {
                    vertex_type.emplace_back(i->point_types[j][k]);
                }
            }
        }
        for(auto i=0;i<critical_points_in_span.size();++i)
        {
            for(auto j=0;j<critical_points_in_span[i].size();++j)
            {
                vertex_type.emplace_back(4);
            }
        }
        

        std::vector<int> edge_type(edges.size(),0);
        for(auto i=0;i<edges.size();++i)
        {            
            
            edge_type[i]=vertex_type[edges[i][0]];
            
        }

        std::ofstream outFile2(edge_value_filename);
        if (!outFile2.is_open()) {
            std::cerr << "Error: Could not open file " << edge_value_filename << " for writing." << std::endl;
            return;
        }
        for (const auto& value : edge_type) {
            outFile2 << value << "\n";  // Each value on a new line
        }
        outFile2.close();
        std::cout << "Edge values saved to " << edge_value_filename << std::endl;
    }


    template<typename T>
    void writeOBJ(const std::string& filename, 
                const std::vector<TraceInSpan<T>>& ridge_valley_in_span,
                const std::vector<std::vector<VectorX<T>>>& critical_points_in_span,
                const std::vector<std::vector<T>>& critical_point_value,
                const std::vector<std::array<size_t,2>>& edges, bool z_zero=false) {


        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }
        // Write vertex data
        if(z_zero)
        {
            for (auto i=ridge_valley_in_span.begin();i<ridge_valley_in_span.end();++i)
            {
                for(auto j=0;j<i->traces.size();++j)
                {
                    for(auto k=0;k<i->traces[j].size();++k)
                    {
                        outFile << std::setprecision(15) << "v " << i->traces[j][k].data()[0] << " " << i->traces[j][k].data()[1] << " " << 0.0 << "\n";
                    }
                }
            }
            for(auto i=0;i<critical_points_in_span.size();++i)
            {
                for(auto j=0;j<critical_points_in_span[i].size();++j)
                {
                    outFile << std::setprecision(15) << "v " << critical_points_in_span[i][j].data()[0] << " " << critical_points_in_span[i][j].data()[1] << " " << 0.0 << "\n";
                }
            }

        }
        else{

            for (auto i=ridge_valley_in_span.begin();i<ridge_valley_in_span.end();++i)
            {
                for(auto j=0;j<i->traces.size();++j)
                {
                    for(auto k=0;k<i->traces[j].size();++k)
                    {
                        outFile<< std::setprecision(15)  << "v " << i->traces[j][k].data()[0] << " " << i->traces[j][k].data()[1] << " " << i->traces_values[j][k] << "\n";
                    }
                }
            }
            for(auto i=0;i<critical_points_in_span.size();++i)
            {
                for(auto j=0;j<critical_points_in_span[i].size();++j)
                {
                    outFile<< std::setprecision(15)  << "v " << critical_points_in_span[i][j].data()[0] << " " << critical_points_in_span[i][j].data()[1] << " " << critical_point_value[i][j] << "\n";
                }
            }
        }


        // Write edge data (lines)
        for (const auto& edge : edges) {
            // OBJ uses 1-based indexing, so add 1 to the indices
            outFile << "l " << edge[0] + 1 << " " << edge[1] + 1 << "\n";
        }
        outFile.close();
        std::cout << "OBJ file written successfully to " << filename << std::endl;

    
    }




    template<typename T>
    bool test_correct(std::vector<TraceInSpan<T>>& trace_in_span)
    {
        for(auto k=0;k<trace_in_span.size();++k)
        {           
            for(auto i=0;i<trace_in_span[k].traces.size();++i)
            {
                for(auto j=0;j<trace_in_span[k].traces[i].size();++j)
                {
                    if(trace_in_span[k].traces[i][j].size()!=2)
                    {
                        std::cout<<"error createing2 "<<trace_in_span[k].traces[i][j].size()<<" "<<j<<std::endl;
                        return false;
                    }
                }
            }
        }

        return true;
    }


    template<typename T>
    void connect_ridge_valley_graph(std::vector<TraceInSpan<T>>& ridge_valley_in_span,
    std::vector<VectorX<T>>& critical_points,
    const Block<T>*b ,size_t real_span_num,
    VectorXi& span_num,
    std::vector<size_t>& valid_span_index, T threshold_square, const std::string& filename,
    T step_size_square, bool z_zero=false, const std::string& edge_value_filename="")
    {
        std::cout<<"threshold_connect_traces_square:"<<threshold_square<<std::endl;
        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);

        //std::vector<size_t>& valid_span_index reverse

        std::unordered_map<size_t, size_t> real_to_valid_span_index;
        real_to_valid_span_index.reserve(real_span_num);

        for (size_t valid_idx = 0; valid_idx < valid_span_index.size(); ++valid_idx) {
            size_t real_idx = valid_span_index[valid_idx];
            real_to_valid_span_index[real_idx] = valid_idx;
        }

       
        std::cout<<"generate real_to_valid_span_index done"<<std::endl;

        std::vector<std::vector<VectorX<T>>> critical_points_in_span;
        std::vector<std::vector<T>> critical_point_value;
        
        if(!critical_points.empty())
        {       
            const VectorX<T> domain_range = b->core_maxs - b->core_mins;
            reorganize_points(b->mfa->var(0).tmesh.all_knots, critical_points, critical_points_in_span, number_in_every_domain,ridge_valley_in_span.size(), real_span_num, b->mfa->var(0).p,real_to_valid_span_index,critical_point_value, b->core_mins, domain_range,valid_span_index);
        }
        for(int i=0;i<ridge_valley_in_span.size();i++)
        {
            ridge_valley_in_span[i].connect_info.resize(ridge_valley_in_span[i].traces.size(),{-1, -1, -1, -1,-1, -1, -1, -1});
        }


        std::cout<<"reorganize critical points in span done"<<std::endl;


        tbb::affinity_partitioner ap;

        int record_valid_span_index = 0;
        int record_valid_span_index_2 = 0;
        std::vector<std::vector<size_t>> around_span(ridge_valley_in_span.size());
        // for(int i=0;i<ridge_valley_in_span.size();i++)
        tbb::parallel_for(tbb::blocked_range<size_t>(0,ridge_valley_in_span.size()),[&](const tbb::blocked_range<size_t>& r)
        {
            for(auto i=r.begin();i<r.end();++i)
            {
                around_span[i].reserve(9);
                find_around_span(valid_span_index[i],number_in_every_domain,span_num,around_span[i],real_to_valid_span_index);
            
            }
        },ap);


        // if(!test_correct(ridge_valley_in_span))
        // {
        //     std::cout<<"error +++6"<<std::endl;
        // }

        // std::cout<<"find around span done"<<std::endl;

        //remove duplicate trace from loops
        for(int i=0;i<ridge_valley_in_span.size();i++)
        {
            delete_duplicate_edges_and_nodes(ridge_valley_in_span[i],around_span[i],step_size_square,ridge_valley_in_span);
        }
        // tbb::parallel_for(tbb::blocked_range<size_t>(0,ridge_valley_in_span.size()),[&](const tbb::blocked_range<size_t>& r)
        // {
        //     for(auto i=r.begin();i<r.end();++i)
        //     {
        //         delete_duplicate_edges_and_nodes(ridge_valley_in_span[i],around_span[i],step_size_square,ridge_valley_in_span);
        //     }
        // },ap);


        std::cout<<"delete duplicate edges and nodes done"<<std::endl;


        size_t prefix_sum=0;
        for(auto i=ridge_valley_in_span.begin();i<ridge_valley_in_span.end();i++)
        {
            i->prefix_start.resize(i->traces.size());
            for(int j=0;j<i->traces.size();j++)
            {
                i->prefix_start[j]=prefix_sum;
                prefix_sum+=i->traces[j].size();
            }             
        }

        //add prefix end of critical points
        std::vector<size_t> critical_points_prefix(critical_points_in_span.size());

        for(auto i=0;i<critical_points_in_span.size();i++)
        {
            if(critical_points_in_span[i].empty())
            {
                continue;
            }
            critical_points_prefix[i]=prefix_sum;
            prefix_sum+=critical_points_in_span[i].size();
        }
        
        T thresholdAngle = 10.0 * M_PI / 180.0;
        T cosThreshold = std::cos(thresholdAngle);

        // can use multi-thread
        T connect_critical_point_threshold =threshold_square;
        if(!critical_points.empty())
        {  
            // for(int i=0;i<ridge_valley_in_span.size();i++)
             tbb::parallel_for(tbb::blocked_range<size_t>(0,ridge_valley_in_span.size()),[&](const tbb::blocked_range<size_t>& r)
            {
                for(auto i=r.begin();i<r.end();++i)
                {
                    compare_with_neighbor_span_critical_point(critical_points_in_span,ridge_valley_in_span[i],i,ridge_valley_in_span,number_in_every_domain,span_num,connect_critical_point_threshold,real_to_valid_span_index,valid_span_index,
                    around_span[i],critical_points_prefix,cosThreshold);
                }
            },ap);

            // for(int i=0;i<ridge_valley_in_span.size();i++)
            // {
            //      compare_with_neighbor_span_critical_point(critical_points_in_span,ridge_valley_in_span[i],i,ridge_valley_in_span,number_in_every_domain,span_num,connect_critical_point_threshold,real_to_valid_span_index,valid_span_index,mfa_data,domain_min,domain_range,
            //         around_span[i],critical_points_prefix,cosThreshold);
            // }

        }

        
        std::cout<<"connect critical points done"<<std::endl;




        
        // // use multi-thread
        // for(int i=0;i<ridge_valley_in_span.size();i++)
        // {
        //     compare_with_neighbor_span_traces(ridge_valley_in_span[i],i,ridge_valley_in_span,threshold_square,around_span[i],step_size_square,critical_points_in_span,cosThreshold);
        // }

        tbb::parallel_for(tbb::blocked_range<size_t>(0,ridge_valley_in_span.size()),[&](const tbb::blocked_range<size_t>& r)
        {
            for(auto i=r.begin();i<r.end();++i)
            {
                compare_with_neighbor_span_traces(ridge_valley_in_span[i],i,ridge_valley_in_span,threshold_square,around_span[i],step_size_square,critical_points_in_span,cosThreshold);
            }
        },ap);
    

        std::cout<<"connect other traces done"<<std::endl;


        // check trace itself.

        tbb::parallel_for(tbb::blocked_range<size_t>(0,ridge_valley_in_span.size()),[&](const tbb::blocked_range<size_t>& r)
        {
            for(auto i=r.begin();i<r.end();++i)
            {
                connect_trace_two_ends(ridge_valley_in_span[i],step_size_square,i);
            }
        },ap);
        
        
        std::cout<<"connect trace two ends done"<<std::endl;

        // for(int i=0;i<ridge_valley_in_span.size();i++)
        // {
        //     complete_connect_info(ridge_valley_in_span[i].connect_info);
        // }

        // std::cout<<ridge_valley_in_span[record_valid_span_index].connect_info[0][4]<<" "<<ridge_valley_in_span[record_valid_span_index].connect_info[0][5]<<" "<<ridge_valley_in_span[record_valid_span_index].connect_info[0][6]<<std::endl;

        T threshold_square_for_adding = threshold_square; //larger to connect point with critical points



        
        // //check point if not connected, connect with critical points
        // if(!critical_points.empty())
        // {           
        //     for(int i=0;i<ridge_valley_in_span.size();i++)
        //     {
        //         compare_with_closed_critical_point(critical_points_in_span,ridge_valley_in_span[i],around_span[i],threshold_square_for_adding);
        //     }
        // }

        // for(int i=0;i<ridge_valley_in_span.size();i++)
        // {
        //     connect_closest_traces(around_span[i],ridge_valley_in_span[i],ridge_valley_in_span,i,threshold_square_for_adding);
        // }

        // for(int i=0;i<ridge_valley_in_span.size();i++)
        // {
        //     complete_connect_info(ridge_valley_in_span[i].connect_info);
        // }

        // std::cout<<ridge_valley_in_span[record_valid_span_index].connect_info[0][4]<<" "<<ridge_valley_in_span[record_valid_span_index].connect_info[0][5]<<" "<<ridge_valley_in_span[record_valid_span_index].connect_info[0][6]<<std::endl;

        //retrieve all edges

        //add prefix start of span


        // std::cout<<ridge_valley_in_span[record_valid_span_index].connect_info[0][4]<<" "<<ridge_valley_in_span[record_valid_span_index].connect_info[0][5]<<" "<<ridge_valley_in_span[record_valid_span_index].connect_info[0][6]<<std::endl;

        // std::cout<<ridge_valley_in_span[record_valid_span_index_2].connect_info[0][0]<<" "<<ridge_valley_in_span[record_valid_span_index_2].connect_info[0][1]<<" "<<ridge_valley_in_span[record_valid_span_index_2].connect_info[0][2]<<std::endl;

        // get edge
        std::vector<std::array<size_t,2>> edge;
        edge.reserve((size_t)(prefix_sum*1.2));



        construct_edge(edge, ridge_valley_in_span,critical_points_in_span,critical_points_prefix,prefix_sum,record_valid_span_index,record_valid_span_index_2);


        std::cout<<"edge size "<<edge.size()<<std::endl;

        writeOBJ(filename,ridge_valley_in_span,critical_points_in_span,critical_point_value,edge,z_zero);
        if(!edge_value_filename.empty())
        {
            writeEdgeValue(edge_value_filename,ridge_valley_in_span,critical_points_in_span,edge);
        }

    }

    //only_read_location is true when only read 2D location
    template<typename T>
    void read_critical_point(std::string& filename, std::vector<VectorX<T>>& critical_points, bool only_read_location=false, bool only_for_a_value=false, T value=0.0, T threshold=1e-2)
    {
        std::vector<Eigen::MatrixXd> root_ori;
        std::vector<Eigen::MatrixXd> root(1);

        std::ifstream file(filename.c_str());
        if (!file) {
            std::cerr << "File does not exist: "<< std::endl;
            return;
        }



        utility::loadMatrixVector(filename.c_str(),root_ori);

        std::cout<<"read critical points 0"<<root_ori[0].rows()<<std::endl;
        

        if(only_read_location)
        {
            root[0]=root_ori[0].block(0,0,root_ori[0].rows(),2);
        }
        else
        {
            root[0]=root_ori[0].block(0,0,root_ori[0].rows(),3);
        }
            

        if(only_for_a_value)
        {
            for(auto i=0;i<root[0].rows();++i)
            {
                if(std::abs(root[0](i,2)-value)<threshold)
                {
                    critical_points.emplace_back(root[0].row(i).transpose());
                }
            }

            std::cout<<"actual critical point size "<<critical_points.size()<<std::endl;
        }
        else
        {
            int size = critical_points.size();
            critical_points.resize(size+root[0].rows());
            for(auto i=0;i<root[0].rows();++i)
            {
                critical_points[i+size]=root[0].row(i).transpose();
            }
        }

    }
  


}