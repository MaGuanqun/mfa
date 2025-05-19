#include <iostream>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <functional>
#include <Eigen/Dense>
#include"parameters.h"

// using namespace Eigen;

namespace spatial_hashing {

template<typename T>
struct VectorHash
{
    size_t operator()(const Eigen::VectorX<T>& v) const
    {
        std::hash<T> hasher;
        size_t seed = 0;
        int index=(v.size()>>1)+1;
        for (int i =1;i<index;++i) {
            seed ^= hasher(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
    
};



// template<typename T>
// struct VectorHashFloor {
//     size_t operator()(const Eigen::VectorX<T>& v) const {
//         std::hash<T> hasher;
//         size_t seed = 0;
//         for (int i =1;i<v.size();++i) {
//             //from boost::hash_combine
//             seed ^= hasher(std::floor(v[i] / v[0]*0.05) )
//             + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//         }
//         return seed;
//     }
// };


// template<typename T>
// struct VectorHashRound {
//     size_t operator()(const Eigen::VectorX<T>& v) const {
//         std::hash<T> hasher;
//         size_t seed = 0;
//         for (int i =1;i<v.size();++i) {
//             //from boost::hash_combine
//             seed ^= hasher(std::round(v[i] / v[0]*0.05) )
//             + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//         }
//         return seed;
//     }
// };


template<typename T>
struct VectorEqual {
    bool operator()(const Eigen::VectorX<T>& a, const Eigen::VectorX<T>& b) const {
        int size=a.size()>>1;
        T squaredDistance = (a.tail(size)-b.tail(size)).squaredNorm();
        return squaredDistance <= a[0] * a[0];
    }
};

template<typename T>
void registerCell(std::vector<Eigen::VectorX<T>>& unique_root, Eigen::VectorX<T>& Position,T epsilon, std::unordered_set<Eigen::VectorX<T>, VectorHash<T>, VectorEqual<T>>& points_step_1,
std::vector<std::vector<size_t>>& temp_index) 
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

    if(Position.size()==5) //dimension = 2
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                Position[1]=i;
                Position[2]=j;

                if(points_step_1.find(Position) != points_step_1.end())
                {
                    return;
                    // deplicated root
                }

                // auto status = points_step_1.insert(Position);
                // if(!status.second)
                // {
                //     duplicated=true;
                // }
            }
        }
    }
    else if(Position.size()==7) //dimension = 3
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                for(auto k:temp_index[2])
                {
                    Position[1]=i;
                    Position[2]=j;
                    Position[3]=k;

                    if(points_step_1.find(Position) != points_step_1.end())
                    {
                        return;
                        // deplicated root
                    }

                    // auto status =  points_step_1.insert(Position);
                    // if(!status.second)
                    // {
                    //     duplicated=true;

                    // }
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
                Position[j+1]=temp_index[j][temp/domain_store_size[j]];
                temp=temp%domain_store_size[j];
            }

            if(points_step_1.find(Position) != points_step_1.end())
            {
                return;
                // deplicated root
            }

            // auto status = points_step_1.insert(Position);
            // if (!status.second)
            // {
            //     duplicated = true;
            // }
        }
    }
    
    // if(!duplicated)
    // {
    unique_root.emplace_back(Position.tail(temp_index.size()));

    // register this point to cells
    if(Position.size()==5) //dimension = 2
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                Position[1]=i;
                Position[2]=j;
                points_step_1.insert(Position);
            }
        }
    }
    else if(Position.size()==7) //dimension = 3
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                for(auto k:temp_index[2])
                {
                    Position[1]=i;
                    Position[2]=j;
                    Position[3]=k;
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
                Position[j+1]=temp_index[j][temp/domain_store_size[j]];
                temp=temp%domain_store_size[j];
            }

            points_step_1.insert(Position);
        }
    }

    // }

}



//use spatial hashing
template<typename T>
void find_all_unique_root(std::vector<Eigen::VectorX<T>>& root, 
std::vector<Eigen::VectorX<T>>& unique_root, T epsilon = SAME_ROOT_EPSILON)
{

    unique_root.reserve(root.size());
    unique_root.clear();
    std::unordered_set<Eigen::VectorX<T>, VectorHash<T>, VectorEqual<T>> points_step_1;

    // std::unordered_set<Eigen::VectorX<T>, VectorHashFloor<T>, VectorEqual<T>> points_step_2;


    points_step_1.reserve(root.size());
    // points_step_2.reserve(root.size());

    Eigen::VectorX<T> temp(2*root[0].size()+1); //[threshold, current hash index, original position]
    temp[0]=epsilon;

    std::vector<std::vector<size_t>> temp_index;
    temp_index.resize(root[0].size());
    for(auto i=temp_index.begin();i<temp_index.end();++i)
    {
        i->reserve(2);
    }

    for(auto i=root.begin();i<root.end();++i)
    {        
        // std::cout<<"*i "<<i->transpose()<<std::endl;
        temp.tail(root[0].size()) = *i;

        registerCell(unique_root,temp,epsilon,points_step_1,temp_index);
    }
    
}

}