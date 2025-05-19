//-----------------
//compute the control points of derivative

//-----------------
#pragma once

#include    <mfa/mfa.hpp>
#include    <mfa/block_base.hpp>


#include    <stdio.h>

#include    <Eigen/Dense>

#include    <random>
#include <tbb/tbb.h>

using namespace std;

using Index = MatrixXf::Index;


// //from write_vtk.cpp PrepSciCtrlPts()
// //compute position of control points
// void ObtainCtrlPts(
//         int&                        nvars,
//         vector< vector <vec3d> >& vars_ctrl_pts,
//         Block<real_t>*              block)
// {
//     vars_ctrl_pts.resize(nvars);
//     vec3d p;

//     for (size_t i = 0; i < nvars; i++)                                      // science variables
//     {
//         // typing shortcuts
//         auto& mfa_data          = block->vars[i].mfa_data;
//         auto& tmesh             = mfa_data->tmesh;
//         auto& tensor_prods      = tmesh.tensor_prods;
//         auto& all_knots         = tmesh.all_knots;
//         auto& all_knot_levels   = tmesh.all_knot_levels;
//         int ndom_dims           = block->geometry.mfa_data->tmesh.tensor_prods[0].ctrl_pts.cols();  // number of geometry dims

//         size_t nctrl_pts = 0;
//         for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
//         {
//             size_t prod = 1;
//             for (auto k = 0; k < ndom_dims; k++)
//                 prod *= tensor_prods[t].nctrl_pts(k);
//             nctrl_pts += prod;
//         }

//         // compute vectors of individual control point coordinates for the tensor product
//         vector<vector<float>> ctrl_pts_coords(ndom_dims);
//         for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
//         {
//             auto& tc = block->vars[i].mfa_data->tmesh.tensor_prods[t];

//             for (auto k = 0; k < ndom_dims; k++)                            // domain dimensions
//             {
//                 for (auto j = 0; j < tc.nctrl_pts(k); j++)                  // control points
//                 {
//                     // offset the knot to the correct control point
//                     KnotIdx knot_min, idx;
//                     if (tc.knot_mins[k] == 0)
//                         knot_min = (mfa_data->p(k) + 1) / 2;
//                     else
//                         knot_min = tc.knot_mins[k];
//                     if (!tmesh.knot_idx_ofst(tc, knot_min, j, k, false, idx))
//                         throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));

//                     float tsum = all_knots[k][idx];                         // odd degree, tsum is on the knot

//                     // odd degree, second control point from global edge is an average
//                     if ((mfa_data->p(k) % 2 == 1 && tc.knot_mins[k] == 0 && j == 1) ||
//                         (mfa_data->p(k) % 2 == 1 && tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2))
//                     {
//                         KnotIdx idx1;
//                         if (tc.knot_mins[k] == 0 && j == 1)
//                         {
//                             if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
//                                 throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));
//                         }
//                         else if (tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2)
//                         {
//                             if (!tmesh.knot_idx_ofst(tc, idx, -1, k, false, idx1))
//                                 throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));
//                         }
//                         tsum += all_knots[k][idx1];
//                         tsum /= 2.0;
//                     }

//                     if (mfa_data->p(k) % 2 == 0)                            // even degree, find center of knot span
//                     {
//                         KnotIdx idx1;
//                         if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
//                             throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));
//                         tsum += all_knots[k][idx1];
//                         tsum /= 2.0;
//                     }
//                     ctrl_pts_coords[k].push_back(block->core_mins(k) + tsum * (block->core_maxs(k) - block->core_mins(k)));

//                     // debug
//                     // fmt::print(stderr, "t {} k {} j {} tsum (param) {} ctrl_pts_coord {}\n",
//                     //         t, k, j, tsum, ctrl_pts_coords[k].back());



//                 }   // control points

//             }   // domain dimensions
//         }   // tensor products

//         // form the tensor product of control points from the vectors of individual coordinates
//         VectorXi ofst = VectorXi::Zero(3);                                              // offset of indices for current tensor
//         for (auto t = 0; t < tmesh.tensor_prods.size(); t++)                            // tensor products
//         {
//             auto& tc   = tmesh.tensor_prods[t];
//             mfa::VolIterator vol_iter(tc.nctrl_pts);
//             VectorXi ijk(ndom_dims);
//             while (!vol_iter.done())                                                    // control points
//             {
//                 vol_iter.idx_ijk(vol_iter.cur_iter(), ijk);

//                 if (tc.weights(vol_iter.cur_iter()) == MFA_NAW)
//                 {
//                     vol_iter.incr_iter();
//                     continue;
//                 }

//                 // first 3 dims stored as mesh geometry
//                 // control point position and optionally science variable, if the total fits in 3d
//                 p.x = ctrl_pts_coords[0][ofst(0) + ijk(0)];
//                 if (ndom_dims < 2)
//                 {
//                     p.y = tc.ctrl_pts(vol_iter.cur_iter(), 0);
//                     p.z = 0.0;
//                 }
//                 else
//                 {
//                     p.y = ctrl_pts_coords[1][ofst(1) + ijk(1)];
//                     if (ndom_dims < 3)
//                         p.z = tc.ctrl_pts(vol_iter.cur_iter(), 0);
//                     else
//                         p.z = ctrl_pts_coords[2][ofst(2) + ijk(2)];
//                 }
//                 vars_ctrl_pts[i].push_back(p);

//                 // science variable also stored as data

//                 // debug
// //                 fmt::print(stderr, "t {} ctrl_pt [{} {} {}]\n",
// //                         t, vars_ctrl_pts[i].back().x, vars_ctrl_pts[i].back().y, vars_ctrl_data[i][vars_ctrl_pts[i].size() - 1]);

//                 vol_iter.incr_iter();
//             }   // control points
//             ofst.head(ndom_dims) += tc.nctrl_pts;
//         }   // tensor products


       

//     }   // science variables
// }






//P98 algorithm 3.3
//curve derivative control points
//find some way avoid copying for control_points
template <typename T>
void CurveDerivCpts(
    const MatrixX<T>&             control_points,
    std::vector<MatrixX<T>>& deriv_control_points,
    const int                     p,//polynominal degree
    const vector<T>&              all_knots,//U
    const int              deriv,
    const T domain_range)
{
    int r = control_points.rows()-1;// =n

    deriv_control_points.resize(deriv + 1);
    //original control points


    // for(int i=0;i<all_knots.size();++i){
    //     std::cout<<all_knots[i]<<" ";
    // }
    // std::cout<<std::endl;

    deriv_control_points[0]=control_points;
    
    //deriv control points
    for(int k=1;k<=deriv;++k)
    {
        deriv_control_points[k].resize(control_points.rows()-k,control_points.cols());
        int tmp=p-k+1;
        for(int i=0;i<=r-k;++i)
        {
            deriv_control_points[k].row(i) = (tmp/(domain_range * (all_knots[i+p+1] - all_knots[i+k])))*(deriv_control_points[k-1].row(i+1) - deriv_control_points[k-1].row(i));

        }


    }


}



template <typename T>
void partial_deriv_control_points(int domain_id, const T* ori_control_points,T* deriv_control_points,
const VectorXi& control_points_number,
const int deriv, const int p,const std::vector<T> &all_knots,T domain_range)
{   
    size_t cpts_number = control_points_number.prod();
    size_t block_size= 1;
    for(int i=0;i<domain_id;++i)
    {
        block_size*=control_points_number[i];
    }
    size_t large_block_size = block_size*control_points_number[domain_id];



    size_t new_ctrl_point_number = control_points_number[domain_id]-deriv;

    auto itr_num=cpts_number/large_block_size;

    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,itr_num),
    [&](const tbb::blocked_range<size_t>& range)
    {
        for(auto m=range.begin();m!=range.end();++m)
        {
            auto i=large_block_size*m;
            size_t new_i = i/control_points_number[domain_id]*new_ctrl_point_number;

            MatrixX<T> cpts;
            cpts.resize(control_points_number[domain_id],1);
            std::vector<MatrixX<T>> temp_deriv_control_points;
            
            for(auto j=0;j<block_size;j++)
            {

                for(auto k=0;k<control_points_number[domain_id];k++)
                {
                    // std::cout<<i<<" "<<j<<" "<<k<<" "<<block_size<<std::endl;
                    cpts.data()[k]=ori_control_points[i+j+k*block_size];
                }
                // std::cout<<"===="<<std::endl;
                // std::cout<<cpts.transpose()<<std::endl;

                CurveDerivCpts(cpts,temp_deriv_control_points,p,all_knots,deriv,domain_range);
                
                // std::cout<<temp_deriv_control_points[deriv].transpose()<<std::endl;

                for(auto k=0;k<new_ctrl_point_number;k++)
                {
                    // std::cout<<new_i<<" "<<j<<" "<<k<<" "<<block_size<<std::endl;
                    deriv_control_points[new_i+j+k*block_size] = temp_deriv_control_points[deriv].data()[k];
                }          
            }

        }
    },ap
    );


}





//(...,d,d,0,0,0) based on (...,d,0,0,0,0)
//I should directly build function (1,0) and control point to compare with the existed result
template <typename T>
void DerivCptsCertainDimension(int domain_id, T* pre_deriv_control_points,T* deriv_control_points,
const VectorXi& control_points_number,
const int deriv, const int p,const std::vector<T> &all_knots,T domain_range)
{   
    VectorXi current_point_number = control_points_number;
    for(int i=0;i<domain_id;++i) 
    {
        current_point_number[i]-=deriv;
    }
    size_t cpts_number = current_point_number.prod();
    size_t block_size= 1;
    for(int i=0;i<domain_id;++i)
    {
        block_size*=current_point_number[i];
    }
    size_t large_block_size = block_size*current_point_number[domain_id];

    MatrixX<T> cpts;
    cpts.resize(current_point_number[domain_id],1);
    std::vector<MatrixX<T>> temp_deriv_control_points;
    size_t new_i;
    size_t new_ctrl_point_number = current_point_number[domain_id]-deriv;
    for(auto i=0;i<cpts_number;i+= large_block_size)
    {
        new_i = i/current_point_number[domain_id]*new_ctrl_point_number;
        for(auto j=0;j<block_size;j++)
        {
            for(auto k=0;k<current_point_number[domain_id];k++)
            {
                // std::cout<<i<<" "<<j<<" "<<k<<" "<<block_size<<std::endl;
                cpts.data()[k]=pre_deriv_control_points[i+j+k*block_size];
            }
            // std::cout<<"===="<<std::endl;
            // std::cout<<cpts.transpose()<<std::endl;

            CurveDerivCpts(cpts,temp_deriv_control_points,p,all_knots,deriv,domain_range);
            
            // std::cout<<temp_deriv_control_points[deriv].transpose()<<std::endl;

            for(auto k=0;k<new_ctrl_point_number;k++)
            {
                // std::cout<<new_i<<" "<<j<<" "<<k<<" "<<block_size<<std::endl;
                deriv_control_points[new_i+j+k*block_size] = temp_deriv_control_points[deriv].data()[k];
            }          
        }
    }

}




//P114 algorithm 3.7
template <typename T>
void DerivCpts(const VectorXi& control_points_number, const std::vector<std::vector<T>>& all_knots,
const VectorXi& p, std::vector<MatrixX<T>>& deriv_control_points, const MatrixX<T>& ctrl_pts, const int deriv,
Eigen::VectorX<T>& domain_range)
{

    // size_t deriv_number=pow((deriv+1),control_points_number.size());
    deriv_control_points.resize(p.size());


    VectorXi deriv_index(p.size());
    deriv_index.setZero();    
    for(auto i=0;i<p.size();++i)
    {
        deriv_index.setZero();    
        deriv_index[i]=deriv;
        size_t size = (control_points_number-deriv_index).prod();
        deriv_control_points[i].resize(size,ctrl_pts.cols());
    }


    size_t cpts_number = control_points_number.prod();


    

    VectorXi derivative_control_point_number = control_points_number;

    //initialization, derivative order (d,0,0)
    derivative_control_point_number[0]=control_points_number[0]-deriv;

    auto itr_num=cpts_number/control_points_number[0];

    tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,itr_num),
    [&](const tbb::blocked_range<size_t>& range)
    {
        for(auto k=range.begin();k!=range.end();++k)
        {

            auto i=control_points_number[0]*k;
            MatrixX<T> ori_control_points;
            ori_control_points = ctrl_pts.block(i,0,control_points_number[0],ctrl_pts.cols());
            std::vector<MatrixX<T>> temp_deriv_control_points;
            CurveDerivCpts(ori_control_points,temp_deriv_control_points,p[0],all_knots[0],deriv,domain_range[0]);
            deriv_control_points[0].block(i/control_points_number[0]*derivative_control_point_number[0],0,derivative_control_point_number[0],ctrl_pts.cols())=temp_deriv_control_points[deriv];
            // std::cout<<ori_control_points.transpose()<<std::endl;
            // std::cout<<temp_deriv_control_points[deriv].transpose()<<std::endl;
            // std::cout<<"==="<<std::endl;

        }
    },ap
    );
    //
    // std::cout<<deriv_control_points[0].transpose()<<std::endl;
    for(int i=1;i<p.size();++i)
    {
        //here assume pt_dim-dom_dim=1
        partial_deriv_control_points(i,ctrl_pts.data(),deriv_control_points[i].data(),
        control_points_number,deriv,p[i],all_knots[i],domain_range[i]);
    }

    // std::cout<<deriv_control_points.size()<<std::endl;
    // std::cout<<deriv_control_points[1]<<std::endl;


}


//compute the x-axis position
template<typename T>
void CurveGeoDerivCpts(
    int&                        nvars,
    Block<T>*              block,
    const VectorXi&         derivs,
    std::vector<std::vector<std::vector<T>>>& ctrl_pts_coords) //[nars,]
{
    ctrl_pts_coords.resize(nvars);
    // for(auto i=control_point.begin();i<control_point.end();++i){
    //     i->resize();
    // }
    for (size_t i = 0; i < nvars; i++)                                      // science variables
    {
        // typing shortcuts
        // auto& mfa_data          = block->vars[i].mfa_data;
        auto& mfa_data          = block->mfa->var(i);
        auto& tmesh             = mfa_data.tmesh;
        auto& tensor_prods      = tmesh.tensor_prods;
        auto& all_knots         = tmesh.all_knots;
        auto& all_knot_levels   = tmesh.all_knot_levels;
        int ndom_dims           = block->dom_dim;  // number of geometry dims

        // size_t nctrl_pts = 0;
        // for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
        // {
        //     size_t prod = 1;
        //     for (auto k = 0; k < ndom_dims; k++)
        //         prod *= (tensor_prods[t].nctrl_pts(k)-derivs[k]);
        //     nctrl_pts += prod;
        // }

        ctrl_pts_coords[i].resize(ndom_dims);
        // compute vectors of individual control point coordinates for the tensor product

        for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
        {
            auto& tc = tmesh.tensor_prods[t];

            for (auto k = 0; k < ndom_dims; k++)                            // domain dimensions
            {

                size_t cpt_size= tc.nctrl_pts(k)-derivs[k];

                ctrl_pts_coords[i][k].resize(cpt_size);


        tbb::affinity_partitioner ap;
        tbb::parallel_for(tbb::blocked_range<size_t>(0,cpt_size),
        [&](const tbb::blocked_range<size_t>& range)
        {
            for(auto j=range.begin();j!=range.end();++j)
            {
                // offset the knot to the correct control point
                KnotIdx knot_min, idx;
                if (tc.knot_mins[k] == 0)
                    knot_min = (mfa_data.p(k)+ 1) / 2;
                else
                    knot_min = tc.knot_mins[k];
                if (!tmesh.knot_idx_ofst(tc, knot_min, j, k, false, idx))
                    throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));
                
                idx+=derivs[k];//for derivative 

                float tsum = all_knots[k][idx];                         // odd degree, tsum is on the knot
                // odd degree, second control point from global edge is an average
                if (((mfa_data.p(k)-derivs[k]) % 2 == 1 && tc.knot_mins[k] == 0 && j == 1) ||
                    ((mfa_data.p(k)-derivs[k]) % 2 == 1 && tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2 - derivs[k]))
                {
                    KnotIdx idx1;
                    if (tc.knot_mins[k] == 0 && j == 1)
                    {
                        if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                            throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));
                    }
                    else if (tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2- derivs[k])
                    {
                        if (!tmesh.knot_idx_ofst(tc, idx, -1, k, false, idx1))
                            throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));
                    }
                    tsum += all_knots[k][idx1];
                    tsum /= 2.0;
                }

                if ((mfa_data.p(k) - derivs[k]) % 2 == 0)                            // even degree, find center of knot span
                {
                    KnotIdx idx1;
                    if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                        throw mfa::MFAError(fmt::format("PrepCtrlPts(): unable to offset knot"));
                    tsum += all_knots[k][idx1];
                    tsum /= 2.0;
                }
                ctrl_pts_coords[i][k][j]=block->core_mins(k) + tsum * (block->core_maxs(k) - block->core_mins(k));
                // debug
                // fmt::print(stderr, "t {} k {} j {} tsum (param) {} ctrl_pts_coord {}\n",
                //         t, k, j, tsum, ctrl_pts_coords[k].back());


            }

        },ap);

            }   // domain dimensions

        
            // for(int m=0;m<ctrl_pts_coords.size();++m)
            // {
            //     std::cout<<ctrl_pts_coords[m].size()<<std::endl;
            // }


        }   // tensor products


        // form the tensor product of control points from the vectors of individual coordinates
        // VectorXi ofst = VectorXi::Zero(3);                                              // offset of indices for current tensor
        // for (auto t = 0; t < tmesh.tensor_prods.size(); t++)                            // tensor products
        // {
        //     auto& tc   = tmesh.tensor_prods[t];
        //     mfa::VolIterator vol_iter(tc.nctrl_pts-derivs);
        //     VectorXi ijk(ndom_dims);
        //     while (!vol_iter.done())                                                    // control points
        //     {
        //         vol_iter.idx_ijk(vol_iter.cur_iter(), ijk);

        //         if (tc.weights(vol_iter.cur_iter()) == MFA_NAW)
        //         {
        //             vol_iter.incr_iter();
        //             continue;
        //         }

        //         for(int k=0;k<ndom_dims;k++)
        //         {
        //             control_point[i][k].push_back(ctrl_pts_coords[k][ofst(k) + ijk(k)]);
        //         }               

        //         vol_iter.incr_iter();
        //     }   // control points
        //     ofst.head(ndom_dims) += tc.nctrl_pts;
        //     ofst.head(ndom_dims) -=derivs;
        // }   // tensor products
       

    }   // science variables


    // for(auto i=control_point[0].begin();i<control_point[0].end();++i)
    // {
    //     std::cout<<*i<<" ";
    // }
    // std::cout<<std::endl;


}


//differentiate of one sci variable in a block
//P114, algorithm 3.7
template <typename T>
void diff_control_point_variable(
    const mfa::MFA_Data<T>& mfa_data,
    int         deriv,
    std::vector<MatrixX<T>>& deriv_control_points,
    Eigen::VectorX<T>& domain_range)
{

    // for(int i=0;i<derivs.size();++i){
    //     CurveDerivCpts(mfa_data.tmesh.tensor_prods[0].ctrl_pts, deriv_control_points,
    //     mfa_data.p[i],mfa_data.tmesh.all_knots[i],derivs[i],domain_range[i]);
    // }

    DerivCpts(mfa_data.tmesh.tensor_prods[0].nctrl_pts,mfa_data.tmesh.all_knots,mfa_data.p,deriv_control_points,
    mfa_data.tmesh.tensor_prods[0].ctrl_pts,deriv,domain_range);    

}


// P114 algorithm 3.7
// surface derivative control points

//Have not finished, even with different parameters, only for first derivative
template <typename T>
void diff_control_point_block(
            Block<real_t>* b, 
            int                               deriv,    // which derivative to take (1 = 1st, 2 = 2nd, ...) in each domain dim.
            std::vector<std::vector<std::vector<std::vector<T>>>>& geo_control_point, //[deriv][vars][dom][...]
            std::vector<std::vector<MatrixX<T>>>& sci_deriv_control_points //[vas][partial_deriv][...]
            )
            // int                               var,
            // int                               partial)  // limit to partial derivative in just this dimension (-1 = no limit))
{

    // vector< vector <vec3d> > vars_ctrl_pts;

    int ndom_dims   = b->dom_dim;          // number of geometry dims
    int nvars           = b->mfa->nvars();                       // number of science variables
    int pt_dim          = b->pt_dim;    

    geo_control_point.resize(deriv+1);

    VectorXi derivs;
    derivs.resize(b->dom_dim);


    for(int i=0;i<=deriv;i++)
    {
        derivs = i * VectorXi::Ones(b->dom_dim);
        CurveGeoDerivCpts(nvars,b,derivs,geo_control_point[i]);
    } 
    
    Eigen::VectorX<T> local_domain_range=b->core_maxs-b->core_mins;

    // for(int i=0;i<geo_control_point[0].size();++i)
    // {
    //     for(int j=0;j<geo_control_point[0][i].size();++j)
    //     {
    //         std::cout<<geo_control_point[0][i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }  
    // std::cin.get();

    // b->approx = new mfa::PointSet<T>(b->input->params, b->input->pt_dim);  // Set decode params from input params

    sci_deriv_control_points.resize(b->mfa->nvars());

    for (auto i = 0; i < b->mfa->nvars(); i++)
        // if (var < 0 || var == i)
        {           
            diff_control_point_variable(b->mfa->var(i), deriv,sci_deriv_control_points[i],local_domain_range); // assumes each variable is scalar
        }
    // combine the control points
    // int domain_dim= b->geometry.mfa_data->tmesh.tensor_prods[0].ctrl_pts.cols();


    // for(int i=0;i<geo_control_point.size();i++)
    // {
    //     for(int j=0;j<geo_control_point[i][0].size();j++)
    //     {
    //         for(int k=0;k<geo_control_point[i][0][j].size();k++)
    //         {   
    //             std::cout<<geo_control_point[i][0][j][k]<<" ";

    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    //     for(int j=0;j<sci_deriv_control_points[0].size();j++)
    //     {
    //         for(int k=0;k<sci_deriv_control_points[0][j].size();k++)
    //         {
    //             std::cout<<sci_deriv_control_points[0][j].data()[k]<<" ";

    //         }
    //         std::cout<<std::endl;
    //     }
    
    // deriv_control_points.resize(b->vars.size());

    // for(int i=0;i<deriv_control_points.size();++i){
    //     deriv_control_points[i].resize(sci_deriv_control_points[i].size());
    //     for(int j=0;j<deriv_control_points[i].size();++j)
    //     {

    //     }



    //     deriv_control_points[i].resize(geo_control_point[i][0].size(),domain_dim+1);

    //     for(int j=0;j<domain_dim;++j){
    //         memcpy(deriv_control_points[i].data()+j*geo_control_point[i][j].size(),geo_control_point[i][j].data(),
    //         sizeof(T)*geo_control_point[i][j].size());
    //     }
    //     memcpy(deriv_control_points[i].data() + domain_dim*geo_control_point[i][0].size(),sci_deriv_control_points[i].data(),
    //         sizeof(T)*geo_control_point[i][0].size());
    // }


}