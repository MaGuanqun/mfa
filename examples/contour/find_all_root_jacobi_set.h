// This is the method using newton method to find all roots
//Hubbard, John, Dierk Schleicher, and Scott Sutherland. "How to find all roots of complex polynomials by Newton’s method." Inventiones mathematicae 146.1 (2001): 1-33.
//https://math.stackexchange.com/questions/998333/finding-the-all-roots-of-a-polynomial-by-using-newton-raphson-method
//This is the version that can only find root in [0,1]
//Only with real number inital points

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"

#include "utility_function.h"

#include "kdtree.h"
#include "mfa_extend.h"

// using namespace std;

#include"parameters.h"


namespace find_all_roots_jacobi_set
{

// Function to find the statistical mode of a list
template<typename T>
int statistical_mode(std::vector<T>& list) {
    if(list.empty()){
        return 1;
    } 
    std::vector<std::pair<int, T>> c;

    std::sort(list.begin(), list.end());
    int length = 1;
    for (size_t i = 1; i < list.size(); ++i) {
        if (list[i] == list[i - 1]) {
            length++;
        } else {
            c.emplace_back(std::make_pair(length, list[i-1]));
            length=1;
        }
    }
    c.emplace_back(std::make_pair(length,list.back()));

    if (c.size() == 1) {
        return c[0].second;
    }

    int mc = 0;
    for (const auto& p : c) {
        if(mc<p.first){
            mc=p.first;
        }       
    }
    for (const auto& p : c) {
        if(mc==p.first){
            return p.second;
        }
    }

    return 1;
    
}




//check all dimensions except the chosen one
template<typename T>
bool InBlockOtherDim(std::vector<std::vector<T>>& span_range, VectorX<T>& point, int dim)
{
    
    for(int i=0;i<span_range.size();++i)
    {
        if(i!=dim)
        {
            if(point[i]<span_range[i][0] || point[i]>span_range[i][1])
            {
                return false;
            }
        }
    }
    return true;
}

template<typename T>
bool findIntersection(std::vector<std::vector<T>>& span_range, VectorX<T>& pre_point, VectorX<T>& current_point,
VectorX<T>& intersection, T root_finding_epsilon)
{
    T t;
    for(int i=0;i<pre_point.size();++i)
    {

        if(abs(current_point[i]-pre_point[i])<root_finding_epsilon)
        {
            continue;
        }

        t = (span_range[i][0]-pre_point[i])/(current_point[i]-pre_point[i]);
        if(t<0.0 || t>1.0)
        {
            t = (span_range[i][1]-pre_point[i])/(current_point[i]-pre_point[i]);

            intersection[i]=span_range[i][1];
        }
        else
        {
            intersection[i]=span_range[i][0];
        }

        for(int j=0;j<pre_point.size();++j)
        {
            if(j!=i)
            {
                intersection[j]=pre_point[j] + t*(current_point[j]-pre_point[j]);
            }
        }

        
        if(InBlockOtherDim(span_range,intersection,i))
        {
            if(t<root_finding_epsilon)
            {
                return false;
            }
            return true;
        }
    }
    
    
    return false;
}



template<typename T>
void compute_h(const mfa::MFA<T>* mfa, const Block<T>* b, VectorX<T>& p0, VectorX<T>& h,
        const mfa::MFA<T>* mfa2, const Block<T>* b2)
{


    VectorX<T> f_first_deriv;
    f_first_deriv.resize(p0.size());
    VectorX<T> dev_f_vector(1);  
    VectorXi deriv(p0.size());
    h.resize(2);
    for(int i=0;i<p0.size();i++)
    {
        deriv.setZero();
        deriv[i]+=1;
        mfa_extend::recover_mfa(mfa,b,p0,dev_f_vector, deriv);
        f_first_deriv[i] = dev_f_vector[0];///local_domain_range[i];
    }
    


    VectorX<T> g_first_deriv;
    g_first_deriv.resize(p0.size());
    for(int i=0;i<p0.size();i++)
    {
        deriv.setZero();
        deriv[i]+=1;
        mfa_extend::recover_mfa(mfa2,b2,p0,dev_f_vector,deriv);
        g_first_deriv[i] = dev_f_vector[0];///local_domain_range[i];
    }

    h[0]=f_first_deriv[0]*g_first_deriv[1]-f_first_deriv[1]*g_first_deriv[0];

}



template<typename T>
void compute_h_dev_h(const mfa::MFA<T>* mfa, const Block<T>* b, VectorXi& span_index, VectorX<T>& p0, VectorX<T>& h, MatrixX<T>& h_first_deriv,
                //MatrixX<T>&             ctrl_pts,   //control points of first derivative
                VectorX<T>&             weights, const mfa::MFA<T>* mfa2, const Block<T>* b2)
{
    VectorX<T> f_first_deriv;
    MatrixX<T> f_second_deriv;
    std::array<T,4> f_third_deriv; 
    //[f_xxx, f_xxy, f_xyy, f_yyy]

    f_second_deriv.resize(span_index.size(),span_index.size());
    f_first_deriv.resize(span_index.size());

    VectorX<T> f_vector(1);
    VectorX<T> dev_f_vector(1);    

    // std::cout<<local_domain_range.transpose()<<std::endl;

    VectorXi deriv(span_index.size());

    for(int i=0;i<span_index.size();i++)
    {
        for(int j=i;j<span_index.size();j++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[j]+=1;
            mfa_extend::recover_mfa(mfa,b,p0,dev_f_vector,deriv);
            f_second_deriv(j,i) = dev_f_vector[0];// / (local_domain_range[i]*local_domain_range[j]);
            f_second_deriv(i,j) = dev_f_vector[0];
        }
    }

    for(int i=0;i<span_index.size();i++)
    {
        deriv.setZero();
        deriv[i]+=1;
        mfa_extend::recover_mfa(mfa,b, p0,f_vector,deriv);
        f_first_deriv[i] = f_vector[0];///local_domain_range[i];
    }


    deriv[0]=3;
    deriv[1]=0;
    for(int i=0;i<4;i++)
    {
        mfa_extend::recover_mfa(mfa,b,p0,dev_f_vector,deriv);
        f_third_deriv[i] = dev_f_vector[0];///local_domain_range[i];
        deriv[0]-=1;
        deriv[1]+=1;
    }

    
    //for g
    VectorX<T> g_first_deriv;
    MatrixX<T> g_second_deriv;
    std::array<T,4> g_third_deriv;

    g_second_deriv.resize(span_index.size(),span_index.size());
    g_first_deriv.resize(span_index.size());

    for(int i=0;i<span_index.size();i++)
    {
        for(int j=i;j<span_index.size();j++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[j]+=1;
            mfa_extend::recover_mfa(mfa2,b2, p0,dev_f_vector,deriv);
            g_second_deriv(j,i) = dev_f_vector[0];// / (local_domain_range[i]*local_domain_range[j]);g
            g_second_deriv(i,j) = dev_f_vector[0];
        }
    }

    for(int i=0;i<span_index.size();i++)
    {
        deriv.setZero();
        deriv[i]+=1;
        mfa_extend::recover_mfa(mfa2,b2, p0,f_vector,deriv);
        g_first_deriv[i] = f_vector[0];///local_domain_range[i];
    }

    deriv[0]=3;
    deriv[1]=0;
    for(int i=0;i<4;i++)
    {
        mfa_extend::recover_mfa(mfa2,b2, p0,dev_f_vector,deriv);
        g_third_deriv[i] = dev_f_vector[0];///local_domain_range[i];
        deriv[0]-=1;
        deriv[1]+=1;
    }

    h_first_deriv.resize(span_index.size(),span_index.size());
    h.resize(span_index.size());
       
    h[0]=f_first_deriv[0]*g_second_deriv.data()[1]+g_first_deriv[1]*f_second_deriv.data()[0]-f_first_deriv[1]*g_second_deriv.data()[0]-g_first_deriv[0]*f_second_deriv.data()[1];
    h[1]=f_first_deriv[0]*g_second_deriv.data()[3]+g_first_deriv[1]*f_second_deriv.data()[1]-f_first_deriv[1]*g_second_deriv.data()[1]-g_first_deriv[0]*f_second_deriv.data()[3];

    h_first_deriv(0,0)=f_second_deriv(0,0)*g_second_deriv(0,1)+f_first_deriv[0]*g_third_deriv[1]+g_second_deriv(0,1)*f_second_deriv(0,0) + g_first_deriv[1]*f_third_deriv[0]-f_second_deriv(0,1)*g_second_deriv(0,0)-f_first_deriv[1]*g_third_deriv[0]- g_second_deriv(0,0)*f_second_deriv(0,1)-g_first_deriv[0]*f_third_deriv[1];

    h_first_deriv(1,0) = h_first_deriv(0,1)=f_second_deriv(0,1)*g_second_deriv(0,1)+f_first_deriv[0]*g_third_deriv[2]+g_second_deriv(1,1)*f_second_deriv(0,0) + g_first_deriv[1]*f_third_deriv[1]-f_second_deriv(1,1)*g_second_deriv(0,0)-f_first_deriv[1]*g_third_deriv[1]- g_second_deriv(1,0)*f_second_deriv(0,1)-g_first_deriv[0]*f_third_deriv[2];

    h_first_deriv(1,1)=f_second_deriv(1,0)*g_second_deriv(1,1)+f_first_deriv[0]*g_third_deriv[3]+g_second_deriv(1,1)*f_second_deriv(1,0) + g_first_deriv[1]*f_third_deriv[2]-f_second_deriv(1,1)*g_second_deriv(1,0)-f_first_deriv[1]*g_third_deriv[2]- g_second_deriv(1,0)*f_second_deriv(1,1)-g_first_deriv[0]*f_third_deriv[3];

}


template<typename T>
void compute_f_dev_f(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data, VectorX<T>& p, VectorX<T>& f, MatrixX<T>& dev_f, Eigen::VectorX<T>& local_domain_range
                )
{

    dev_f.resize(p.size(),p.size());
    f.resize(p.size());

    VectorX<T> f_vector(p.size());
    VectorX<T> dev_f_vector(p.size());    
    // std::cout<<local_domain_range.transpose()<<std::endl;

    VectorXi deriv(p.size());

    for(int i=0;i<p.size();i++)
    {
        for(int j=i;j<p.size();j++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[j]+=1;
            mfa->DecodePt(*mfa_data,p,deriv,dev_f_vector);
            dev_f(j,i) =dev_f_vector[0];// / (local_domain_range[i]*local_domain_range[j]);
            dev_f(i,j) = dev_f_vector[0];
        
        }
    }

    for(int i=0;i<p.size();i++)
    {
        deriv.setZero();
        deriv[i]+=1;
        mfa->DecodePt(*mfa_data,p,deriv,f_vector);
        f[i] = f_vector[0];///local_domain_range[i];
    }

}



//only use itr as convergence condition
template<typename T>
void newton_itr(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data, std::vector<int>& span, std::vector<T>& result,T p, int max_itr,
                //MatrixX<T>&             ctrl_pts,   //control points of first derivative
                VectorX<T>&             weights,Eigen::VectorXd& local_domain_range)
{
    result.clear();
    T f,dev_f;
    int i=0;
    while (i<max_itr)
    {
        compute_f_dev_f(mfa,mfa_data, span, p,f,dev_f,weights,local_domain_range);
        p=p-f/dev_f;
        result.emplace_back(p);
        i++;
    }
    
}


// newton method with single initial_point
template<typename T>
bool newton(const mfa::MFA<T>* mfa, const Block<T>* b,VectorXi& span_index, VectorX<T>& result, VectorX<T>& p, int max_itr,
               // MatrixX<T>&             ctrl_pts,   //control points of first derivative
                VectorX<T>&             weights,
                std::vector<std::vector<T>>& span_range,
                T d_max_square, VectorX<T>& center,
                bool& filtered_out,Eigen::VectorX<T>& local_domain_range,Eigen::VectorX<T>& local_min, size_t& itr_num, T root_finding_epsilon, const mfa::MFA<T>* mfa2, const Block<T>* b2, Eigen::VectorX<T>& domain_min,Eigen::VectorX<T>& domain_range)
{
    itr_num=0;


    MatrixX<T> dev_h;
    VectorX<T> h;
    compute_h_dev_h(mfa,b,span_index,p,h,dev_h,weights, mfa2, b2);




    if(h.squaredNorm()<root_finding_epsilon*root_finding_epsilon)
    {
        result = p;
        return true;
    }

    result=p;

    T temp_rec;

    VectorX<T> pre_point = p;
    VectorX<T> intersection = p;

    while(itr_num<max_itr)
    {


        if(h.size()==2)
        {
            double determinant = dev_h.data()[3]*dev_h.data()[0]-dev_h.data()[1]*dev_h.data()[2];
            if(std::abs(determinant) < HESSIAN_DET_EPSILON_JACOBI_SET)
            {
                return false;                
            }
            else
            {  
                MatrixX<T> inv(2,2);
                inv.data()[0]=dev_h.data()[3];
                inv.data()[1]=-dev_h.data()[1];
                inv.data()[2]=-dev_h.data()[2];
                inv.data()[3]=dev_h.data()[0];           
                inv/=determinant;
                p -= inv*h;
            }
          
        }
        else
        {
            // printf("determinant: %f\n",dev_f.determinant());
            // std::cout<<p.transpose()<<std::endl;
            // std::cout<<"p=="<<std::endl;
            // std::cout<<f.transpose()<<std::endl;
            // std::cout<<"f=="<<std::endl;
            
            if(std::abs(dev_h.determinant()) < HESSIAN_DET_EPSILON_JACOBI_SET)
            {
                return false;
            }

            p -= dev_h.colPivHouseholderQr().solve(h);                  
        }


        // std::cout<<"p "<<p.transpose()<<std::endl;


        // if(initial_p[0]>0.5&&initial_p[0]<0.51
        // && initial_p[1]>0.5 && initial_p[1]<0.51)
        // {
        //     std::cout<<p.transpose()<<std::endl;
        // }

        if((p-center).squaredNorm()>d_max_square)
        {
            filtered_out = true;
            return false;
        }
        if(!utility::InDomain(p,local_min,local_domain_range))
        {
            filtered_out = true;
            return false;
        }

        // if(!InBlock(span_range,p))
        // {
        //     // if(findIntersection(span_range,pre_point,p,intersection))
        //     // {
        //     //     p=intersection;
        //     //     pre_point = intersection;
        //     // }    
        //     return false;
        // }


        compute_h_dev_h(mfa,b,span_index,p,h,dev_h,weights, mfa2, b2);   


        if(itr_num>0){
            if(h.squaredNorm()<root_finding_epsilon*root_finding_epsilon){//|| (result-p).squaredNorm()<ROOT_FINDING_EPSILON*ROOT_FINDING_EPSILON
              
                if(!utility::InBlock(span_range,p))
                {
                    return false;
                }
                
                               

                result = p;
                return true;
            }
        }

        // result = p;
        itr_num++;
    }

    return false;

}




template<typename T>
bool newRoot(VectorX<T>& z, std::vector<VectorX<T>>& root_so_far, T threshold)
{
    for(int i=0;i<root_so_far.size();++i)
    {
        if((z-root_so_far[i]).squaredNorm()<threshold*threshold)
        {
            return false;
        }
    }
    return true;
}


// Function to find the roots of the polynomial using Newton's method
template<typename T>
void newtonSolve(const mfa::MFA<T>* mfa, const mfa::MFA_Data<T>& mfa_data, const Block<T>* b, VectorXi& span_index, std::vector<VectorX<T>>& root,
    VectorX<T>&             weights, int& original_root_size, 
    size_t& filtered_out_num,Eigen::VectorX<T>& local_domain_range, Eigen::VectorX<T>& local_min, T same_root_threshold, size_t& itr_num, T root_finding_epsilon, int maxIter, const mfa::MFA<T>* mfa2, const mfa::MFA_Data<T>& mfa_data2, const Block<T>* b2, 
    VectorX<T>& domain_min, VectorX<T>& domain_range) { 

    VectorXi one = VectorXi::Ones( mfa_data.p.size());
    // int deg = (mfa_data->p-one).prod();


    

    std::vector<VectorX<T>> root_record;
    
    // std::cout<<"max_iteration--"<<maxIter<<std::endl;

    std::vector<std::vector<T>> span_range(span_index.size());


    VectorX<T> center(span_index.size());
    for(int i=0;i<span_index.size();++i)
    {    
        span_range[i].emplace_back(mfa_data.tmesh.all_knots[i][span_index[i]]*local_domain_range[i]+local_min[i]);
        span_range[i].emplace_back(mfa_data.tmesh.all_knots[i][span_index[i]+1]*local_domain_range[i]+local_min[i]);

        center[i]=(span_range[i][0]+span_range[i][1])*0.5;
    }   



    // compute distance to terminate iteration
    T d_max_square=0;
    for(auto i=span_range.begin();i!=span_range.end();++i)
    {  
        d_max_square+=((*i)[1]-(*i)[0])*((*i)[1]-(*i)[0]);
    }

        
    d_max_square*=DISTRANCE_STOP_ITR * DISTRANCE_STOP_ITR; //d^2=(2*diagonal of span)^2

    std::vector<std::vector<T>>initial_point;


    // VectorX<T> p; p.resize(3);
    // p[0] =84.53500366210938;
    // p[1] =76.00250244140625;
    // p[2] = 88.05249786376953;
    // VectorX<T> ini_p = (p-local_min).cwiseQuotient(local_domain_range);
    // std::vector<T> ini_p_vector;
    // ini_p_vector.emplace_back(ini_p[0]);
    // ini_p_vector.emplace_back(ini_p[1]);
    // ini_p_vector.emplace_back(ini_p[2]);
    // initial_point.emplace_back(ini_p_vector);
    // std::cout<<p.transpose()<<std::endl;
    // std::cout<<"ini_p "<<ini_p.transpose()<<std::endl;


    utility::compute_initial_points_js(initial_point,mfa_data.p,span_range);

    VectorXi num_initial_point_every_domain(initial_point.size());
    for(int i=0;i<num_initial_point_every_domain.size();i++)
    {
        num_initial_point_every_domain[i]=initial_point[i].size();
    }

    int num_initial_point = num_initial_point_every_domain.prod();

    VectorX<T> next_root; 

    int total_root_num=0;




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


    bool filetered_out = false;
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
        filetered_out = false;

        size_t current_itr_num=0;
        bool find_root=newton(mfa, b, span_index, next_root, current_initial_point,maxIter,weights,span_range,d_max_square,center,filetered_out,local_domain_range,local_min,current_itr_num,root_finding_epsilon,mfa2, b2,domain_min,domain_range);
        
        itr_num += current_itr_num;

        if(find_root)
        {
            
            // std::cout<<"is a new root "<<std::endl;
     
            
                original_root_size++;
                if(newRoot(next_root,root_record,same_root_threshold))
                {       
                    root.emplace_back(next_root);      
                    root_record.emplace_back(next_root);

                    total_root_num+=1;              

                    // std::cout<<next_root_in_original_domain.transpose()<<std::endl;
                    // std::cout<<"record+++++"<<std::endl;
                    // if(total_root_num>=deg)
                    // {
                    //     break;
                    // }
                } 

            // }

            
        }

        if(filetered_out)
        {
            filtered_out_num++;
        }
        
        // break;
        
    }

}


template<typename T>
void newtonSolve(Block<T>* block, std::vector<std::vector<VectorXi>>& span_index, 
std::vector<VectorX<T>>& root,//std::vector<int>& multi_of_root,
    //MatrixX<T>&             ctrl_pts,   //control points of first derivative
    VectorX<T>&             weights, int current_index,
    int& original_root_size, size_t& filtered_out_num, T same_root_threshold, size_t& itr_num, T root_finding_epsilon, int maxItr,Block<T>* block2) //(p+1)^d initial points) 
{
    VectorX<T> local_domain_range = block->core_maxs-block->core_mins;

    for(auto i=0;i<block->mfa->nvars();++i)
    {
        newtonSolve(block->mfa,block->mfa->var(i), block, span_index[i][current_index], root,weights, 
        original_root_size,
       filtered_out_num,local_domain_range,block->core_mins,same_root_threshold,itr_num,root_finding_epsilon,maxItr,block2->mfa,block2->mfa->var(i),block2,block->core_mins,local_domain_range);
    }
}

// //the original root are in [0,1], convert it to 
// template<typename T>
// void convertToDomain(VectorX<T>& core_min,VectorX<T>& range,std::vector<std::vector<VectorX<T>>>& ori_root,std::vector<std::vector<VectorX<T>>>& domain_root)
// {
//     domain_root.clear();
//     domain_root.resize(ori_root.size());

//     tbb::affinity_partitioner ap;
//     tbb::parallel_for(tbb::blocked_range<size_t>(0,ori_root.size()),
//     [&](const tbb::blocked_range<size_t>& interval)
//     {
//         for(auto i=interval.begin();i!=interval.end();++i)
//         {
//             domain_root[i].resize(ori_root[i].size());
//             for(int j=0;j<ori_root[i].size();++j)
//             {
//                 domain_root[i][j]=core_min+ ori_root[i][j].cwiseProduct(range);

//             }
//         }
//     },ap
//     );

// }


//the original root are in [0,1], convert it to 
template<typename T>
void convertToDomain(VectorX<T>& core_min,VectorX<T>& range,std::vector<std::vector<VectorX<T>>>& ori_root,std::vector<std::vector<VectorX<T>>>& domain_root)
{
    domain_root.clear();
    domain_root.resize(ori_root.size());

    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,ori_root.size()),
    [&](const tbb::blocked_range<size_t>& interval)
    {
        for(auto i=interval.begin();i!=interval.end();++i)
        {
            domain_root[i].resize(ori_root[i].size());
            for(int j=0;j<ori_root[i].size();++j)
            {
                domain_root[i][j]=core_min+ ori_root[i][j].cwiseProduct(range);

            }
        }
    },ap
    );

}

template<typename T>
void convertRootToMatrix(std::vector<std::vector<VectorX<T>>>& root_vec_in_domain, MatrixXd& record_root, std::vector<std::vector<T>>& func_value)
{
    size_t total_root_num=0;
    for(auto i=0;i<root_vec_in_domain.size();++i)
    {
        total_root_num+=root_vec_in_domain[i].size();
    }
    record_root.resize(total_root_num,root_vec_in_domain[0].size()+1);
    size_t index=0;
    for(auto i=0;i<root_vec_in_domain.size();++i)
    {
        for(auto j=0;j<root_vec_in_domain[i].size();++j)
        {
            record_root.block(index+j,0,1,root_vec_in_domain[0][0].size())=root_vec_in_domain[i][j].transpose();
            record_root(index+j,root_vec_in_domain[0][0].size())=func_value[i][j];
        }
        index+=root_vec_in_domain[i].size();
    }
}

template<typename T>
void convertToDomain(VectorX<T>& core_min,VectorX<T>& range,std::vector<VectorX<T>>& ori_root,std::vector<VectorX<T>>& domain_root)
{
    domain_root.clear();
    domain_root.resize(ori_root.size());

    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,ori_root.size()),
    [&](const tbb::blocked_range<size_t>& interval)
    {
        for(auto i=interval.begin();i!=interval.end();++i)
        {
            domain_root[i]=core_min+ ori_root[i].cwiseProduct(range);
        }
    },ap
    );

}



template<typename T>
void convertFromDomain(std::vector<VectorX<T>>& domain_root, std::vector<VectorX<T>>& uniform_root,VectorX<T>& core_min,VectorX<T>& range)
{

    uniform_root.clear();
    uniform_root.resize(domain_root.size());

    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,domain_root.size()),
    [&](const tbb::blocked_range<size_t>& interval)
    {
        for(auto i=interval.begin();i!=interval.end();++i)
        {
            uniform_root[i]=(domain_root[i]-core_min).cwiseQuotient(range);
        }
    },ap
    );
}





//only need critical points with function value 0
template<typename T>
void filterRoot(const mfa::MFA<T>* mfa, const Block<T>* b, std::vector<VectorX<T>>& root,std::vector<VectorX<T>>& new_root, std::vector<VectorX<T>>& value0, T threshold,const mfa::MFA<T>* mfa2, const Block<T>* b2)
{
    std::vector<VectorX<T>> value;
    value.resize(root.size());

    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, root.size()),
        [&](const tbb::blocked_range<size_t>& range)
        {
            for(auto i=range.begin();i!=range.end();++i)
            {

                compute_h(mfa,b,root[i],value[i],mfa2,b);

            }
        },ap
    );

    for(auto i=0;i<root.size();++i)
    {
        if(std::abs(value[i][0])<threshold)
        {
            new_root.emplace_back(root[i]);
            value0.emplace_back(value[i]);
        }
    }

}

//only need critical points with function value 0
template<typename T>
void functionValue0(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data, std::vector<VectorX<T>>& root, std::vector<VectorX<T>>& domain_root, std::vector<VectorX<T>>& new_domain_root, std::vector<VectorX<T>>& new_root, std::vector<VectorX<T>>& value0, T threshold)
{
    std::vector<VectorX<T>> value;
    value.resize(root.size());

    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, root.size()),
        [&](const tbb::blocked_range<size_t>& range)
        {
            for(auto i=range.begin();i!=range.end();++i)
            {
                mfa->DecodePt(*mfa_data,root[i],value[i]);

            }
        },ap
    );

    for(auto i=0;i<root.size();++i)
    {
        new_root.emplace_back(root[i]);
        value0.emplace_back(value[i]);
        new_domain_root.emplace_back(domain_root[i]);
        
    }

}


template<typename T>
int compute_index(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data,VectorX<T>& root, T threshold)
{
    VectorXi deriv(root.size());
    MatrixX<T> Hessian;

    Hessian.resize(root.size(),root.size());

    VectorX<T> deriv_value(root.size());    

    for(int i=0;i<root.size();i++)
    {
        for(int j=i;j<root.size();j++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[j]+=1;
            mfa->DecodePt(*mfa_data,root,deriv,deriv_value); 
            Hessian(j,i) = deriv_value[0];// / (local_domain_range[i]*local_domain_range[j]);
            Hessian(i,j) = deriv_value[0];
        }
    }

    if(abs(Hessian.determinant()) < threshold)
    {
        return -1;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<T>> solver(Hessian);

    // Eigenvalues are guaranteed to be real
    Eigen::VectorX<T> eigenvalues = solver.eigenvalues();
    int index=0;
    for (int i = 0; i < eigenvalues.size(); i++)
    {
        if(eigenvalues[i]<0)
        {
            index++;
        }
    }

    return index;

}



template<typename T>
void get_critical_point_index(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data, std::vector<VectorX<T>>& root,std::vector<int>& index, T threshold)
{
    index.resize(root.size());
    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,index.size()),
    [&](const tbb::blocked_range<size_t>& range)
    {
        for(auto i=range.begin();i!=range.end();++i)
        {
            index[i]=compute_index(mfa,mfa_data,root[i],threshold);
        }

    },ap
    );
}

template<typename T>
void getDerivative(mfa::MFA<T>* mfa, mfa::MFA_Data<T>* mfa_data, std::vector<VectorX<T>>& root,std::vector<VectorX<T>>& value)
{
    value.resize(root.size());

    int size = root[0].size();
    VectorXi deriv(size);
    VectorX<T> deriv_value(size);

    for(auto i =0;i<root.size();++i)
    {
        value[i].resize(size);
        // std::cout<<root[i].transpose()<<std::endl;
        for(int j=0;j<size;j++)
        {
            deriv.setZero();
            deriv[j]=1;
            mfa->DecodePt(*mfa_data,root[i],deriv,deriv_value);   
            value[i][j]=deriv_value[0];        
        }
        // std::cout<<root[i].transpose()<<" "<<value[i][0].transpose()<<" "<<value[i][1].transpose()<<std::endl;
    }

}

template<typename T>
T getAverageNorm(std::vector<VectorX<T>>& deriv)
{
    T sum=0;
    for(auto i=0;i<deriv.size();++i)
    {
        sum+=deriv[i].norm();
    }
    sum/=deriv.size();
    return sum;

}


// template<typename T>
// bool VectorXEqual(const std::vector<T>& vec1, const std::vector<T>& vec2)
// {
//     // Define the equality criterion based on distance
//     return (vec1 - vec2).squaredNorm() < RATIO_CHECK_SAME_ROOT*RATIO_CHECK_SAME_ROOT*ROOT_FINDING_EPSILON*ROOT_FINDING_EPSILON;
// }




// template<typename T>
// bool VectorXSmaller(const Eigen::VectorX<T>& vec1, const std::vector<T>& vec2)
// {
//     for(int i=vec2.size()-1;i>=0;i--)
//     {
//         if(vec1[i]<vec2[i])
//         {
//             return true;
//         }
//         if(vec1[i]>vec2[i])
//         {
//             return false;
//         }

//     }
//     return false;
// }


// find all root1 elements in root2 
template<typename T>
void find_all_overlapped_root(std::vector<VectorX<T>>& root1,std::vector<VectorX<T>>& root2, std::vector<VectorX<T>>& overlapped_root,
T accuracy)
{
    overlapped_root.clear();
    KDTree kd(root2);
    std::vector<int> neighbor_points;
    for(int i=0;i<root1.size();i++)
    {
        neighbor_points.clear();
        kd.radiusSearch(root2, root1[i],neighbor_points,accuracy*accuracy);
        if(!neighbor_points.empty())
        {
            overlapped_root.emplace_back(root1[i]);
        }
    }

}



//use kdtree
template<typename T>
void find_all_unique_root(std::vector<VectorX<T>>& root, std::vector<VectorX<T>>& unique_root,
T accuracy, std::vector<int>&duplicated_number)
{
    unique_root.clear();
    KDTree kd(root);
    std::vector<bool> merged(root.size(),false);

    std::vector<int> neighbor_points;
    neighbor_points.reserve(std::pow(2,root[0].size()));

    // duplicated_number.reserve(root.size());

    for(int i=0;i<root.size();i++)
    {
        if(!merged[i])
        {
            neighbor_points.clear();
            kd.radiusSearch(root,root[i],neighbor_points,accuracy*accuracy);
            for(auto j=neighbor_points.begin();j!=neighbor_points.end();j++)
            {
                if((*j)!=i)
                {
                    merged[*j]=true;
                }
            }
            // duplicated_number.emplace_back(neighbor_points.size());
        }
    }

    unique_root.reserve(root.size());

    for(int i=0;i<root.size();++i)
    {
        if(!merged[i])
        {
            unique_root.emplace_back(root[i]);
        }
    }


    unique_root.shrink_to_fit();
    // duplicated_number.shrink_to_fit();


    // if(duplicated_number.size()!=unique_root.size())
    // {
    //     std::cout<<"error: the number of unique root not compatible with duplicated number recording"<<std::endl;
    // }
}

template<typename T>
void create_initial_point_in_a_range(Eigen::VectorX<T> point, std::vector<Eigen::VectorX<T>>& initial_point, T half_cube, int point_num_each_dim)
{

    // std::cout<<half_cube<<std::endl; 

    std::vector<std::vector<T>> initial_points_every_domain(point.size());
    for(int i=0;i<point.size();i++)
    {
        initial_points_every_domain[i].reserve(point_num_each_dim);
        for (int j = 0; j < point_num_each_dim; ++j) {
            T coe = ((T(j))/(point_num_each_dim-1))*(half_cube+half_cube);
            initial_points_every_domain[i].emplace_back(point[i]-half_cube+coe);
        }
    }

    
    VectorXi num_initial_point_every_domain(initial_points_every_domain.size());
    for(int i=0;i<num_initial_point_every_domain.size();i++)
    {
        num_initial_point_every_domain[i]=initial_points_every_domain[i].size();
    }

    int num_initial_point = num_initial_point_every_domain.prod();

    VectorXi domain_index;
    VectorXi number_in_every_domain;
    VectorX<T> current_initial_point(initial_points_every_domain.size());
    utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);

    for(int i=0;i<num_initial_point;++i)
    {
        //std::cout<<"test for wrong "<<std::endl;

        utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
       // std::cout<<i<<" "<<domain_index.transpose()<<std::endl;
        for(int j=0;j<current_initial_point.size();j++)
        {
            current_initial_point[j]=initial_points_every_domain[j][domain_index[j]];
        }        
        initial_point.emplace_back(current_initial_point);
    }

}



    template<typename T>
    void newton_method(Block<real_t>* b,std::vector<std::vector<VectorX<T>>>& different_root_1_from_2, T root_finding_epsilon)
    {
        int control_points_num = b->mfa->var(0).tmesh.tensor_prods[0].nctrl_pts(0);
        std::cout<<" control_points_num "<<control_points_num<<std::endl;
        std::cout<< b->mfa->var(0).p.transpose()<<std::endl;
        T range = 1.0/double(control_points_num - b->mfa->var(0).p[0]);
        int maxIter = 200;

        std::vector<std::vector<T>> span_range(different_root_1_from_2[0][0].size());
        for(int i=0;i<span_range.size();i++)
        {
            span_range[i].emplace_back(0);
            span_range[i].emplace_back(1);
        }
        Eigen::VectorX<T> local_domain_range=b->core_maxs-b->core_mins;

        for(int i=0;i<different_root_1_from_2[0].size();i++)
        {
            std::vector<VectorX<T>> root;
            std::vector<Eigen::VectorX<T>> initial_point;
            VectorX<T> next_root; 
            create_initial_point_in_a_range(different_root_1_from_2[0][i],initial_point,range,9);

            // for(int j=0;j<initial_point.size();j++)
            // {
            //     std::cout<<initial_point[j].transpose()<<std::endl;
            // }
            // std::cout<<"finished creating_initial points "<<std::endl;

            for(int j=0;j<initial_point.size();++j)
            {
                if(newton(b->mfa, b->mfa->var(0),next_root, initial_point[j],maxIter,span_range,local_domain_range,root_finding_epsilon))
                {
                    if(newRoot(next_root,root,SAME_ROOT_EPSILON))
                    {       
                        root.emplace_back(next_root);                      

                    } 
                    
                }
            }

            bool has_connection = false;
            for(int j=0;j<root.size();j++)
            {
                if((root[j]-different_root_1_from_2[0][i]).squaredNorm()<1e-4)
                {
                    has_connection = true;
                    std::cout<<(b->core_mins + root[j].cwiseProduct(b->core_maxs-b->core_mins)).transpose()<<" "<<(b->core_mins + different_root_1_from_2[0][i].cwiseProduct(b->core_maxs-b->core_mins)).transpose()<<std::endl;
                    break;
                }
            }

            if(!has_connection)
            {
                std::cout<<"cannot find a nearby root for "<<different_root_1_from_2[0][i].transpose()<<std::endl;
            }

        }

    }


}