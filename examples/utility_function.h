#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <Eigen/Dense>




namespace utility
{
    template<class T>
    size_t obtain_index_from_domain_index(T& domain_index, VectorXi& number_in_every_domain)//column major number_in_every_domain:[1,p0,p0*p1,p0*p1*p2,...]
    {
        size_t index = domain_index[0];
        for(int k=1;k<number_in_every_domain.size();k++)
        {
            index += number_in_every_domain[k]*domain_index[k];
        }
        return index;
    }

    template<class T>
    size_t obtain_index_from_domain_index_row_major(T& domain_index, VectorXi& number_in_every_domain)//column major number_in_every_domain:[...,p0*p1*p2,p0*p1,p0,1]
    {
        size_t index = domain_index[domain_index.size()-1];
        for(int k=number_in_every_domain.size()-2;k>=0;k--)
        {
            index += number_in_every_domain[k]*domain_index[k];
        }
        return index;
    }


    template<class T>
    void obtain_number_in_every_domain(T& ori_num,VectorXi& number_in_every_domain)
    {
        number_in_every_domain.resize(ori_num.size());
        number_in_every_domain[0]=1;
        for(int k=1;k<ori_num.size();k++)
        {
            number_in_every_domain[k]=number_in_every_domain[k-1]*ori_num[k-1];
        }
    }

    template<class T>
    void obtain_number_in_every_domain_row_major(T& ori_num,VectorXi& number_in_every_domain)
    {
        number_in_every_domain.resize(ori_num.size());
        number_in_every_domain[number_in_every_domain.size()-1]=1;
        for(int k=ori_num.size()-2;k>=0;k--)
        {
            number_in_every_domain[k]=number_in_every_domain[k+1]*ori_num[k+1];
        }
    }


    void obtainDomainIndex(size_t index, VectorXi& domain_index,
    VectorXi& number_in_every_domain) //number_in_every_domain:[1,p0,p0*p1,p0*p1*p2,...]
    {
        domain_index.resize(number_in_every_domain.size());
        size_t temp_0, temp_1;

        size_t current_index = index;
        size_t temp_index=index;
        for(auto i=number_in_every_domain.size()-1;i>0;i--)
        {
            domain_index[i] = temp_index / number_in_every_domain[i];
            temp_index = temp_index % number_in_every_domain[i];
        }
        domain_index[0]=temp_index; 
    }

    template<class T>
    void loadMatrixVector(const char* filename,std::vector<MatrixX<T>>& root)
    {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        int vector_size;
        in.read((char*)(&vector_size), sizeof(int));
        root.resize(vector_size);
        size_t matrix_col, matrix_row;
        for(int i=0;i<vector_size;i++)
        {
            in.read((char*)(&matrix_row), sizeof(size_t));
            in.read((char*)(&matrix_col), sizeof(size_t));
            root[i].resize(matrix_row,matrix_col);
            in.read((char*)(root[i].data()), matrix_row*matrix_col*sizeof(double));
        }
        in.close();
    }

    template<class T>
    void writeMatrixVector(const char* filename,std::vector<MatrixX<T>>& root)
    {
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);        
        int vector_size;
        vector_size = root.size();
        out.write((char*)(&vector_size), sizeof(int));
        size_t matrix_col, matrix_row;
        for(int i=0;i<vector_size;i++)
        {
            matrix_row = root[i].rows();
            matrix_col = root[i].cols();
            out.write((char*)(&matrix_row), sizeof(size_t));
            out.write((char*)(&matrix_col), sizeof(size_t));
            out.write((char*)(root[i].data()), matrix_row*matrix_col*sizeof(T)); 
        }
        out.close();
    }

    template<class T>
    void saveToCSV(const std::string& filename, const std::vector<Eigen::VectorX<T>>& vec, T function_value) {

        std::ofstream file(filename);
        
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        for(int i=0;i<vec[0].size();i++)
        {
            file << "x" << i << ",";
        }
        file << "value\n";
        
        for (const auto& vector : vec) {
            for (int i = 0; i < vector.size(); ++i) {
                file << vector[i];
                file << ",";            
            }
            file << function_value;
            file << "\n";
        }
        
        file.close();
    }


    template<class T>
    void saveToCSV(const std::string& filename, const std::vector<Eigen::VectorX<T>>& vec) {

        std::ofstream file(filename);
        
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        for(int i=0;i<vec[0].size()-1;i++)
        {
            file << "x" << i << ",";
        }
        file << "x"<<vec[0].size()-1<< "\n";
        
        for (const auto& vector : vec) {
            for (int i = 0; i < vector.size() - 1; ++i) {
                file << vector[i];
                file << ",";            
            }
            file << vector[vector.size() - 1];
            file << "\n";
        }
        
        file.close();
    }


    template<class T>
    bool InDomain(VectorX<T>& point) // point should lie in [0,1]^d
    {
        for(int i=0;i<point.size();++i)
        {
            if(point[i]<0 || point[i]>1)
            {
                return false;
            }
        }
        return true;

    }

    template<class T>
    bool InDomain(VectorX<T>& point,const VectorX<T>&         domain_min,
            const VectorX<T>&         domain_range) // point should lie in [0,1]^d
    {
        for(int i=0;i<point.size();++i)
        {
            if(point[i]<domain_min[i] || (point[i] - domain_min[i])>domain_range[i])
            {
                return false;
            }
        }
        return true;

    }

    
    template<typename T>
    bool InBlock(std::vector<std::vector<T>>& span_range, VectorX<T>& point) //span_range: [[min0,max0],[min1,max1]]
    {
        
        for(int i=0;i<span_range.size();++i)
        {
            // if(point[i]<=0 || point[i]>=1)
            if(point[i]<span_range[i][0] || point[i]>span_range[i][1])
            {
                return false;
            }
        }
        return true;
    }

    //checkl if cross the zero
    template <typename T>
    bool check_valid_span(T* value, std::vector<size_t>& index_for_control_point, T func_value)
    {
        T max, min;
        max = value[index_for_control_point[0]];
        min = value[index_for_control_point[0]];
        for(auto i=index_for_control_point.begin()+1;i<index_for_control_point.end();++i)
        {
            if (value[*i] > max)
                max =value[*i];
            if (value[*i] < min)
                min = value[*i];

            
            if(max>=func_value-(1e-5) && min <=func_value+(1e-5))
            {
                return true;
            }
        }

        // std::cout<<max<<" "<<min<<std::endl;

        return false;

    }

    //checkl if cross the zero
    template <typename T>
    bool check_valid_span(T* value, std::vector<size_t>& index_for_control_point)
    {
        T max, min;
        max = value[index_for_control_point[0]];
        min = value[index_for_control_point[0]];
        for(auto i=index_for_control_point.begin()+1;i<index_for_control_point.end();++i)
        {
            if (value[*i] > max)
                max =value[*i];
            if (value[*i] < min)
                min = value[*i];

            
            if(max>=-1e-5 && min <=1e-5)
            {
                return true;
            }
        }

        // std::cout<<max<<" "<<min<<std::endl;

        return false;

    }

    
    template<typename T>
    std::vector<T> compute_initial_points(int d, std::vector<T>& span_range) {
        int s = std::ceil(0.26632 * std::log(d));
        int n = std::ceil(8.32547 * d * std::log(d))*s/2;
        // T r=(1.0 + std::sqrt(2.0))*std::pow(1.0 - 1.0 / d, 1.0 / (4.0 * s));
        std::vector<T> result;
        result.reserve(n);
        for (int j = 0; j < n; ++j) {
            T coe = (T(j)/(n-1))*(span_range[1]-span_range[0]);
            result.emplace_back(span_range[0]+coe);
        }
        return result;
    }


    template<typename T>
    void compute_initial_points(std::vector<T>&initial_points_every_domain, int point_num,  std::vector<T>& span_range) { 
      
        initial_points_every_domain.reserve(point_num);
        for (int j = 0; j < point_num; ++j) {
            T coe = ((0.5+T(j))/(point_num))*(span_range[1]-span_range[0]);
            initial_points_every_domain.emplace_back(span_range[0]+coe);
        }

    }

    //there will be 2^n+1 initial points
    template<typename T>
    void compute_initial_points(std::vector<std::vector<T>>&initial_points_every_domain, const VectorXi& degree,  std::vector<std::vector<T>>& span_range) { 
        // int point_num=std::pow(2,n)+1;
        initial_points_every_domain.resize(span_range.size());
        for(int i=0;i<span_range.size();i++)
        {
            int point_num=degree[i]+1;
            initial_points_every_domain[i].reserve(point_num);
            for (int j = 0; j < point_num; ++j) {
                T coe = ((0.5+T(j))/(point_num))*(span_range[i][1]-span_range[i][0]);
                initial_points_every_domain[i].emplace_back(span_range[i][0]+coe);
            }
            
        }

    }


    //there will be 2^n+1 initial points
    template<typename T>
    void compute_initial_points_iso(std::vector<std::vector<T>>&initial_points_every_domain, VectorXi& degree,  std::vector<std::vector<T>>& span_range) { 
        // int point_num=std::pow(2,n)+1;
        initial_points_every_domain.resize(span_range.size());
        for(int i=0;i<span_range.size();i++)
        {
            int point_num=degree[i]+3;
            initial_points_every_domain[i].reserve(point_num);
            for (int j = 0; j < point_num; ++j) {
                T coe = (T(j)/(point_num-1))*(span_range[i][1]-span_range[i][0]);
                initial_points_every_domain[i].emplace_back(span_range[i][0]+coe);
            }
            
        }

    }

        template<typename T>
    void compute_initial_points_rv(std::vector<std::vector<T>>&initial_points_every_domain, VectorXi& degree,  std::vector<std::vector<T>>& span_range) { 
        // int point_num=std::pow(2,n)+1;
        initial_points_every_domain.resize(span_range.size());
        for(int i=0;i<span_range.size();i++)
        {
            int point_num=(degree[i]*3-1)+3;
            initial_points_every_domain[i].reserve(point_num);
            for (int j = 0; j < point_num; ++j) {
                T coe = (T(j)/(point_num-1)*0.98+0.01)*(span_range[i][1]-span_range[i][0]);
                initial_points_every_domain[i].emplace_back(span_range[i][0]+coe);
            }
            
        }

    }

    
    //jacobi set
        template<typename T> 
    void compute_initial_points_js(std::vector<std::vector<T>>&initial_points_every_domain, VectorXi& degree,  std::vector<std::vector<T>>& span_range) { 
        // int point_num=std::pow(2,n)+1;
        initial_points_every_domain.resize(span_range.size());
        for(int i=0;i<span_range.size();i++)
        {
            int point_num=(degree[i]+degree[i])+2;
            initial_points_every_domain[i].reserve(point_num);
            for (int j = 0; j < point_num; ++j) {
                T coe = (T(j)/(point_num-1)*0.98+0.01)*(span_range[i][1]-span_range[i][0]);
                initial_points_every_domain[i].emplace_back(span_range[i][0]+coe);
            }
            
        }

    }


    void obtain_index_for_control_point(std::vector<size_t>&index_for_control_point,
    VectorXi& ori_span_index, VectorXi& number_in_every_domain_block,
    VectorXi& number_in_every_domain)
    {
        int total_num = index_for_control_point.size();
        VectorXi block_span_index(ori_span_index.size());
        VectorXi current_span_index;
        for(int i=0;i<total_num;++i)
        {
            obtainDomainIndex(i,block_span_index,number_in_every_domain_block);
            current_span_index = ori_span_index+block_span_index;

            index_for_control_point[i]=obtain_index_from_domain_index(current_span_index,number_in_every_domain);

        }

    }

    void obtain_index_for_control_point_row_major(std::vector<size_t>&index_for_control_point,
    VectorXi& ori_span_index, VectorXi& number_in_every_domain_block,
    VectorXi& number_in_every_domain)
    {
        int total_num = index_for_control_point.size();
        VectorXi block_span_index(ori_span_index.size());
        VectorXi current_span_index;
        for(int i=0;i<total_num;++i)
        {
            obtainDomainIndex(i,block_span_index,number_in_every_domain_block);
            current_span_index = ori_span_index+block_span_index;

            index_for_control_point[i]=obtain_index_from_domain_index_row_major(current_span_index,number_in_every_domain);

        }

    }


    template<typename T>
    bool check_if_a_trace_is_in_a_loop(std::vector<VectorX<T>>& loop, std::vector<VectorX<T>>& trace, T threshold_square)
    {
        if(trace.size()>loop.size())
        {
            return false;
        }
        


        for(int i=0;i<loop.size();++i)
        {
            // if(loop[i].size()!=trace[0].size())
            // {
            //     std::cout<<"error "<<loop[i].size()<<" "<<loop.size()<<" "<<trace[0].size()<<std::endl;
            //     for(int k=0;k<loop.size();++k)
            //     {
            //         std::cout<<loop[k].size()<<" ";
            //     }
            // }
            if((loop[i]-trace[0]).squaredNorm()<threshold_square)
            {
                bool find_same_initial = false;
                for(int j=0;j<loop.size();++j)
                {
                    // if(loop[j].size()!=trace.back().size())
                    // {
                    //     std::cout<<"error size "<<loop.size()<<" "<<j<<" " <<loop[j].size()<<" "<<trace.back().size()<<std::endl;
                    //     for(int k=0;k<loop.size();++k)
                    //     {
                    //         std::cout<<loop[k].transpose()<<" ";
                    //     }
                    // }

                   if((loop[j]-trace.back()).squaredNorm()<threshold_square)
                   {
                          find_same_initial = true;
                          break;
                   }
                }

                if(find_same_initial)
                {
                    if(trace.empty())
                    {
                        std::cout<<"empty trace"<<std::endl;
                        return true;
                    }
                    int mid = trace.size()>>1;
                    for(int j=0;j<loop.size();++j)
                    {

                        if((loop[j]-trace[mid]).squaredNorm()<threshold_square)
                        {
                            return true;
                        }
                    }
                }
            }
        }
        
        return false;

    }




}