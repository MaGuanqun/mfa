#include <iostream>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <functional>
#include <Eigen/Dense>
#include <unordered_map>

// using namespace Eigen;

#ifndef FIND_SAME_CENTER_H
#define FIND_SAME_CENTER_H

namespace same_center {

template<typename T>
struct VHash
{
    size_t operator()(const Eigen::VectorX<T>& v) const
    {
        std::hash<T> hasher;
        size_t seed = 0;
        int index=(v.size()>>1)+1;
        for (int i =2;i<index;++i) {
            seed ^= hasher(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
    
};


template<typename T>
struct VEqual {
    bool operator()(const Eigen::VectorX<T>& a, const Eigen::VectorX<T>& b) const {
        int size=(a.size()>>1)-1;
        T squaredDistance = (a.tail(size)-b.tail(size)).squaredNorm();
        return squaredDistance <= a[0] * a[0];
    }
};

template<typename T>
bool findDuplicateCell(Eigen::VectorX<T>& Position,T epsilon, std::unordered_set<Eigen::VectorX<T>, VHash<T>, VEqual<T>>& points_step_1,
std::vector<std::vector<size_t>>& temp_index, size_t& unique_index) 
{
    size_t k1,k2;
    int j=0;
    for(int i=((Position.size()>>1)+1);i<Position.size();++i)
    {
        temp_index[j].clear();
        k1=std::floor(Position[i] / (20*epsilon));
        k2=std::ceil(Position[i] / (20*epsilon));
        temp_index[j].emplace_back(k1); 
        if(k1!=k2)
            temp_index[j].emplace_back(k2);

        j++;
    }

    if(Position.size()==6) //dimension = 2
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                Position[2]=i;
                Position[3]=j;

                auto it=points_step_1.find(Position);
                if( it!= points_step_1.end())
                {
                    unique_index = it->data()[1];
                    return true;
                    // deplicated root
                }
            }
        }
    }
    else if(Position.size()==8) //dimension = 3
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                for(auto k:temp_index[2])
                {
                    Position[2]=i;
                    Position[3]=j;
                    Position[4]=k;
                    auto it=points_step_1.find(Position);
                    if(it != points_step_1.end())
                    {
                        unique_index = it->data()[1];                        
                        return true;
                        // deplicated root
                    }

                }
            }
        }
    }
    else
    {
        int total_case=1;
        for(auto i:temp_index)
        {
            total_case*=i.size();
        }
        std::vector<int> domain_store_size(temp_index.size());
        domain_store_size.back()=1; //domain_store_size [...,d1*d2*d3,d1*d2,d1,1] dn is 1 or 2
        for(int i=temp_index.size()-2;i>=0;--i)
        {
            domain_store_size[i]=domain_store_size[i+1]*temp_index[i+1].size();
        }

        for(int i=0;i<total_case;++i)
        {
            int temp=i;
            for(int j=0;j<temp_index.size();++j)
            {
                Position[j+2]=temp_index[j][temp/domain_store_size[j]];
                temp=temp%domain_store_size[j];
            }
            auto it=points_step_1.find(Position);
            if(it != points_step_1.end())
            {
                unique_index = it->data()[1];
                return true;
                // deplicated root
            }
        }
    }    

    // register this point to cells
    if(Position.size()==6) //dimension = 2
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                Position[2]=i;
                Position[3]=j;
                points_step_1.insert(Position);
            }
        }
    }
    else if(Position.size()==8) //dimension = 3
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                for(auto k:temp_index[2])
                {
                    Position[2]=i;
                    Position[3]=j;
                    Position[4]=k;
                    points_step_1.insert(Position);

                }
            }
        }
    }
    else
    {
        int total_case=1;
        for(auto i:temp_index)
        {
            total_case*=i.size();
        }
        std::vector<int> domain_store_size(temp_index.size());
        domain_store_size.back()=1; //domain_store_size [...,d1*d2*d3,d1*d2,d1,1] dn is 1 or 2
        for(int i=temp_index.size()-2;i>=0;--i)
        {
            domain_store_size[i]=domain_store_size[i+1]*temp_index[i+1].size();
        }

        for(int i=0;i<total_case;++i)
        {
            int temp=i;
            for(int j=0;j<temp_index.size();++j)
            {
                Position[j+2]=temp_index[j][temp/domain_store_size[j]];
                temp=temp%domain_store_size[j];
            }

            points_step_1.insert(Position);
        }
    }

    return false;

}



//use spatial hashing
template<typename T>
void find_all_unique_root(std::vector<Eigen::VectorX<T>>& root, 
std::unordered_map<size_t, std::vector<size_t>>& cluster,
T epsilon)
{
    // points_step_2.reserve(root.size());

    
    // ==========
    std::unordered_set<Eigen::VectorX<T>, VHash<T>, VEqual<T>> points_step_1;
    //[threshold, index_in_unique_vector, current hash index, original position]
    
   
    Eigen::VectorX<T> temp(2*root[0].size()+2);
    temp[0]=epsilon;
    std::vector<std::vector<size_t>> temp_index(root[0].size());
    size_t unique_index;

    for(auto i=0;i<root.size();++i)
    {        
         
         temp[1]=i;
         temp.tail(root[0].size())=root[i];

        if(findDuplicateCell(temp,epsilon,points_step_1,temp_index,unique_index))
        {
            cluster[unique_index].emplace_back(i);
        }
        else
        {
            cluster[i].emplace_back(i);
        }
    }

    // for (const auto& pair : cluster) {
    //     // std::cout << "Key: " << pair.first << " Values: ";
    //     for (int val : pair.second) {

    //         std::cout << val << " ";
    //     }
    //     std::cout << "\n";
    // }
    
}

}


#endif