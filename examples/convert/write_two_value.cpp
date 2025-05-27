//--------------------------------------------------------------
// writes all vtk files for initial, approximated, and control points
//
// optionally generates test data for analytical functions and writes to vtk
//
// output precision is float irrespective whether input is float or double
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#include    "mfa/mfa.hpp"
#include    <iostream>
#include    <stdio.h>

#include    <diy/master.hpp>
#include    <diy/io/block.hpp>

#include    "opts.h"

#include    "writer.hpp"
#include    "block.hpp"

#include <variant>
#include <tbb/tbb.h>
#include "utility_function.h"
#include "write_to_ply.h"
#include "../contour/ridge_valley_graph.h"
#include "mfa_extend.h"


// namespace {
//     tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
// }


template<typename T>
void write_function_pointset_vtk(mfa::PointSet<T>* ps, char* filename,Block<real_t>* block,Block<real_t>* block2,std::chrono::duration<double>& run_time)
{


    if (ps == nullptr)
    {
        cout << "Did not write " << filename << " due to uninitialized pointset" << endl;
        return;
    }
    if (ps->npts == 0)
    {
        cout << "Did not write " << filename << " due to empty pointset" << endl;
        return;
    }

    int dom_dim = ps->dom_dim;
    // int pt_dim  = ps->pt_dim;
    int nvars = 2;   // TODO: this assumes all vars are scalar

    int ori_nvars = 1;


    auto start_time = std::chrono::high_resolution_clock::now();


    vector<int> npts_dim;  // only used if data is structured
    if (ps->is_structured())
    {
        for (size_t k = 0; k < 3; k++)
        {
            if (k < dom_dim) 
                npts_dim.push_back(ps->ndom_pts(k));
            else
                npts_dim.push_back(1);
        }
    }

    float** pt_data = new float*[nvars];
    for (size_t j = 0; j < nvars; j++)
    {
        pt_data[j]  = new float[ps->npts];
    }

    vector<vec3d>   pt_coords;
    pt_coords.resize(ps->npts);

    VectorXd range=block->core_maxs-block->core_mins;


    tbb::affinity_partitioner ap;
    tbb::parallel_for((tbb::blocked_range<size_t>(0,ps->npts)),
    [&](const tbb::blocked_range<size_t>& interval)
    {
        for(size_t j=interval.begin();j<interval.end();++j)
        {
            VectorX<T> result_value(block->mfa->nvars());
            VectorX<T> dom_coordinate = ps->domain.block(j,0,1,dom_dim).transpose();
            //covert geo_coordinate to [0,1]
            // VectorX<T> geo_coordinate = (dom_coordinate-block->core_mins).cwiseQuotient(range);

            for(int m=0;m<dom_dim;m++)
            {
                if(dom_coordinate[m]<block->core_mins[m]){
                    dom_coordinate[m] = block->core_mins[m];
                }
                if(dom_coordinate[m]>block->core_maxs[m]){
                    dom_coordinate[m] = block->core_maxs[m];
                }
            }



            VectorX<T> f_first_deriv(dom_coordinate.size());
            for (int k = 0; k < ori_nvars; k++)                         // science variables
            {
                mfa_extend::recover_mfa(block, dom_coordinate,result_value);
                pt_data[k][j] = result_value[0];
                mfa_extend::recover_mfa(block2, dom_coordinate,result_value);
                pt_data[k+ori_nvars][j] = result_value[0];

            }


            vec3d           pt;
            pt.x = ps->domain(j, 0);
            
            if(dom_dim > 1){
                pt.y = ps->domain(j, 1);
                if(dom_dim<3)
                {
                    pt.z = 0.0; 
                }
                else{
                    pt.z = ps->domain(j, 2);
                }
            }
            else
            {
                pt.y = pt_data[0][j];
                pt.z=0.0;
            }   
            pt_coords[j]=pt;
        }
        
    },ap
    );

    auto end_time = std::chrono::high_resolution_clock::now();

    run_time=end_time-start_time; // measure the time it took to compute the values

    // science variable settings
    int* vardims        = new int[nvars];
    char** varnames     = new char*[nvars];
    int* centerings     = new int[nvars];
    for (int i = 0; i < nvars; i++)
    {
        vardims[i]      = 1;                                // TODO; treating each variable as a scalar (for now)
        varnames[i]     = new char[256];
        centerings[i]   = 1;
        sprintf(varnames[i], "var%d", i);
    }
   // write raw original points
    if (ps->is_structured())
    {
        write_curvilinear_mesh(
            /* const char *filename */                  filename,
            /* int useBinary */                         0,
            /* int *dims */                             &npts_dim[0],
            /* float *pts */                            &(pt_coords[0].x),
            /* int nvars */                             nvars,
            /* int *vardim */                           vardims,
            /* int *centering */                        centerings,
            /* const char * const *varnames */          varnames,
            /* float **vars */                          pt_data);
    }
    else
    {
        write_point_mesh(
        /* const char *filename */                      filename,
        /* int useBinary */                             0,
        /* int npts */                                  pt_coords.size(),
        /* float *pts */                                &(pt_coords[0].x),
        /* int nvars */                                 nvars,
        /* int *vardim */                               vardims,
        /* const char * const *varnames */              varnames,
        /* float **vars */                              pt_data);
    }

    delete[] vardims;
    for (int i = 0; i < nvars; i++)
        delete[] varnames[i];
    delete[] varnames;
    delete[] centerings;
    for (int j = 0; j < nvars; j++)
    {
        delete[] pt_data[j];
    }
    delete[] pt_data;
}



// TODO: Only scalar-valued and 3D vector-valued variables are supported (because of the VTK writer)
// If a variable has a different output dimension, the writer will skip that variable and continue.
template<typename T>
void write_pointset_vtk(mfa::PointSet<T>* ps, char* filename, int sci_var = -1)
{
    if (ps == nullptr)
    {
        cout << "Did not write " << filename << " due to uninitialized pointset" << endl;
        return;
    }
    if (ps->npts == 0)
    {
        cout << "Did not write " << filename << " due to empty pointset" << endl;
        return;
    }

    int dom_dim = ps->dom_dim;
    int geom_dim = ps->geom_dim();
    int nvars = ps->nvars();
    bool include_var = true;        // Include the specified science variable in the geometry coordinates
    int var_col = ps->model_dims().head(sci_var + 1).sum(); // column of the variable to be visualized

    // Sanity checks and modify 'include_var' if settings conflict
    if (geom_dim < 1 || geom_dim > 3)
    {
        cerr << "Did not write " << filename << " due to improper dimension in pointset" << endl;
        return;
    }
    if (sci_var < 0)
    {
        include_var = false;
    }
    else if (ps->var_dim(sci_var) != 1 && geom_dim < 3)
    {
        cerr << "For " << filename << ", specified science variable (#" << sci_var << ") is not a scalar. Output will be planar." << endl;
        include_var = false;
    }

    vector<int> npts_dim;  // only used if data is structured
    if (ps->is_structured())
    {
        for (size_t k = 0; k < 3; k++)
        {
            if (k < dom_dim) 
                npts_dim.push_back(ps->ndom_pts(k));
            else
                npts_dim.push_back(1);
        }
    }

    float** pt_data = new float*[nvars];
    for (size_t k = 0; k < nvars; k++)
    {
        pt_data[k]  = new float[ps->npts * ps->var_dim(k)];
    }

    vec3d           pt;
    vector<vec3d>   pt_coords;
    for (int j = 0; j < ps->npts; j++)
    {
        // Add geometric coordinates
        if (geom_dim == 1)
        {
            pt.x = ps->domain(j, 0);
            pt.y = include_var ? ps->domain(j, var_col) : 0.0;
            pt.z = 0.0;
        }
        else if (geom_dim == 2)
        {
            pt.x = ps->domain(j, 0);
            pt.y = ps->domain(j, 1);
            pt.z = include_var ? ps->domain(j, var_col) : 0.0;
        }
        else
        {
            pt.x = ps->domain(j, 0);
            pt.y = ps->domain(j, 1);
            pt.z = ps->domain(j, 2);
        }
        pt_coords.push_back(pt);

        // Add science variable data
        int offset_idx = 0;
        for (int k = 0; k < nvars; k++)
        {
            int vd = ps->var_dim(k);
            for (int l = 0; l < vd; l++)
            {
                pt_data[k][j*vd + l] = ps->domain(j, geom_dim + offset_idx);
                offset_idx++;
            }
        }    
    }

    // science variable settings
    int* vardims        = new int[nvars];
    char** varnames     = new char*[nvars];
    int* centerings     = new int[nvars];
    for (int k = 0; k < nvars; k++)
    {
        vardims[k]      = ps->var_dim(k);
        varnames[k]     = new char[256];
        centerings[k]   = 1;
        snprintf(varnames[k], 256, "var%d", k);
    }

    // write raw original points
    if (ps->is_structured())
    {
        write_curvilinear_mesh(
            /* const char *filename */                  filename,
            /* int useBinary */                         0,
            /* int *dims */                             &npts_dim[0],
            /* float *pts */                            &(pt_coords[0].x),
            /* int nvars */                             nvars,
            /* int *vardim */                           vardims,
            /* int *centering */                        centerings,
            /* const char * const *varnames */          varnames,
            /* float **vars */                          pt_data);
    }
    else
    {
        write_point_mesh(
        /* const char *filename */                      filename,
        /* int useBinary */                             0,
        /* int npts */                                  pt_coords.size(),
        /* float *pts */                                &(pt_coords[0].x),
        /* int nvars */                                 nvars,
        /* int *vardim */                               vardims,
        /* const char * const *varnames */              varnames,
        /* float **vars */                              pt_data);
    }

    delete[] vardims;
    for (int i = 0; i < nvars; i++)
        delete[] varnames[i];
    delete[] varnames;
    delete[] centerings;
    for (int j = 0; j < nvars; j++)
    {
        delete[] pt_data[j];
    }
    delete[] pt_data;
}


// make combinations of min, max corner vertices in index and real space
void CellVertices(
        int             ndom_dims,                      // number of domain dimensions
        vec3d&          min,                            // min corner
        vec3d&          max,                            // max corner
        vector<vec3d>&  tensor_pts)                     // (output) vertices
{
    vec3d p;

    p.x = min.x;
    p.y = min.y;
    p.z = min.z;
    tensor_pts.push_back(p);

    p.x = max.x;
    tensor_pts.push_back(p);

    if (ndom_dims > 1)
    {
        p.y = max.y;
        tensor_pts.push_back(p);

        p.x = min.x;
        tensor_pts.push_back(p);

        if (ndom_dims > 2)
        {
            p.x = min.x;
            p.y = min.y;
            p.z = max.z;
            tensor_pts.push_back(p);

            p.x = max.x;
            tensor_pts.push_back(p);

            p.y = max.y;
            tensor_pts.push_back(p);

            p.x = min.x;
            tensor_pts.push_back(p);
        }
    }
}

// prep tmesh tensor extents
void PrepTmeshTensorExtents(
    vector<vec3d>&  tensor_pts_real,                    // (output) points in real space
    vector<vec3d>&  tensor_pts_index,                   // (output) points in index space
    vector<int>&    ntensor_pts,
    Block<real_t>*  block)                              // curent block
{
int dom_dim = block->mfa->dom_dim;
int nvars = block->mfa->nvars();

// tmesh tensor extents
ntensor_pts.resize(3);
for (auto i = 0; i < 3; i++)
    ntensor_pts[i] = (i >= dom_dim ? 1 : 2);

for (auto j = 0; j < nvars; j++)
{
    vec3d min_real, max_real, min_index, max_index;         // extents in real and index space

    const mfa::Tmesh<real_t>& tmesh = block->mfa->var(j).tmesh;
    auto& tensor_prods      = tmesh.tensor_prods;
    auto& all_knots         = tmesh.all_knots;
    auto& all_knot_levels   = tmesh.all_knot_levels;

    // form extents in index and real space
    for (auto k = 0; k < tensor_prods.size(); k++)
    {
        min_index.x = tensor_prods[k].knot_mins[0];
        min_real.x  = block->core_mins[0] + all_knots[0][tensor_prods[k].knot_mins[0]] *
            (block->core_maxs[0] - block->core_mins[0]);
        if (dom_dim > 1)
        {
            min_index.y = tensor_prods[k].knot_mins[1];
            min_real.y  = block->core_mins[1] + all_knots[1][tensor_prods[k].knot_mins[1]] *
                (block->core_maxs[1] - block->core_mins[1]);
        }
        else
        {
            min_index.y = 0.0;
            min_real.y  = 0.0;
        }
        if (dom_dim > 2)
        {
            min_index.z = tensor_prods[k].knot_mins[2];
            min_real.z  = block->core_mins[2] + all_knots[2][tensor_prods[k].knot_mins[2]] *
                (block->core_maxs[2] - block->core_mins[2]);
        }
        else
        {
            min_index.z = 0.0;
            min_real.z  = 0.0;
        }

        max_index.x = tensor_prods[k].knot_maxs[0];
        max_real.x  = block->core_mins[0] + all_knots[0][tensor_prods[k].knot_maxs[0]] *
            (block->core_maxs[0] - block->core_mins[0]);
        if (dom_dim > 1)
        {
            max_index.y = tensor_prods[k].knot_maxs[1];
            max_real.y  = block->core_mins[1] + all_knots[1][tensor_prods[k].knot_maxs[1]] *
                (block->core_maxs[1] - block->core_mins[1]);
        }
        else
        {
            max_index.y = 0.0;
            max_real.y  = 0.0;
        }
        if (dom_dim > 2)
        {
            max_index.z = tensor_prods[k].knot_maxs[2];
            max_real.z  = block->core_mins[2] + all_knots[2][tensor_prods[k].knot_maxs[2]] *
                (block->core_maxs[2] - block->core_mins[2]);
        }
        else
        {
            max_index.z = 0.0;
            max_real.z  = 0.0;
        }

        // make vertex points for cells
        CellVertices(dom_dim, min_index, max_index, tensor_pts_index);
        CellVertices(dom_dim, min_real, max_real, tensor_pts_real);

        // debug
//             fmt::print(stderr, "tensor {} extents: index [{} {} {} : {} {} {}] real [{} {} {} : {} {} {}]\n",
//                 k, min_index.x, min_index.y, min_index.z, max_index.x, max_index.y, max_index.z,
//                 min_real.x, min_real.y, min_real.z, max_real.x, max_real.y, max_real.z);

    }   // tensor products
}   // nvars
}

// old version of geometry control points that averages over a window of knots
// based on Youssef's python code
// not used anymore, but kept for reference
void oldPrepGeomCtrlPts(
    vector<vec3d>&              geom_ctrl_pts,
    Block<real_t>*              block)
{
// typing shortcuts
const auto& mfa_data          = block->mfa->geom();
const auto& tmesh             = mfa_data.tmesh;
const auto& tensor_prods      = tmesh.tensor_prods;
const auto& all_knots         = tmesh.all_knots;
const auto& all_knot_levels   = tmesh.all_knot_levels;

int dom_dim     = block->mfa->dom_dim;
int geom_dim    = block->mfa->geom_dim();                                                 // number of geometry dims
vec3d p;


// compute vectors of individual control point coordinates for the tensor product
vector<vector<float>> ctrl_pts_coords(dom_dim);
for (auto t = 0; t < tensor_prods.size(); t++)                                                      // tensor products
{
    const TensorProduct<real_t>& tc = tensor_prods[t];
    for (auto k = 0; k < dom_dim; k++)                                                            // domain dimensions
    {
        int skip = 0;
        // starting knot in sequence for computing control point coordinate
        KnotIdx knot_min = tc.knot_mins[k];
        if (knot_min)
        {
            // skip knots at a deeper level than the tensor
            for (auto l = 0; l < block->mfa->geom().p(k); l++)
            {
                while (all_knot_levels[k][knot_min - l - skip] > tc.level)
                    skip++;
            }
            knot_min -= (block->mfa->geom().p(k) - 1 + skip);
        }

        for (auto j = 0; j < tc.nctrl_pts(k); j++)                      // control points
        {
            float tsum  = 0.0;
            int skip1   = skip;                                         // number of knots at a deeper level that should be skipped
            // skip knots at a deeper level than the tensor
            for (int l = 1; l < block->mfa->geom().p(k) + 1; l++)
            {
                // skip knots at a deeper level than the tensor
                while (all_knot_levels[k][knot_min + j + l + skip1] > tc.level)
                    skip1++;
                tsum += all_knots[k][knot_min + j + l + skip1];
            }
            tsum /= float(block->mfa->geom().p(k));
            ctrl_pts_coords[k].push_back(block->core_mins(k) + tsum * (block->core_maxs(k) - block->core_mins(k)));

            // debug
//                 fprintf(stderr, "t=%d k=%d j=%d tsum=%.3lf ctrl_pts_coord=%.3lf\n", t, k, j, tsum, ctrl_pts_coords[k].back());
        }   // control points
    }   // domain dimensions
}   // tensor products

// form the tensor product of control points from the vectors of individual coordinates
VectorXi ofst = VectorXi::Zero(3);                              // offset of indices for current tensor
for (auto t = 0; t < tensor_prods.size(); t++)                  // tensor products
{
    const TensorProduct<real_t>& tc   = tensor_prods[t];
    mfa::VolIterator vol_iter(tc.nctrl_pts);
    VectorXi ijk(dom_dim);
    while (!vol_iter.done())                                    // control points
    {
        vol_iter.idx_ijk(vol_iter.cur_iter(), ijk);

        if (tc.weights(vol_iter.cur_iter()) == MFA_NAW)
        {
            vol_iter.incr_iter();
            continue;
        }

        // first 3 dims stored as mesh geometry
        p.x = ctrl_pts_coords[0][ofst(0) + ijk(0)];
        if (geom_dim < 2)
            p.y = 0.0;
        else
            p.y = ctrl_pts_coords[1][ofst(1) + ijk(1)];
        if (geom_dim < 3)
            p.z = 0.0;
        else
            p.z = ctrl_pts_coords[2][ofst(2) + ijk(2)];
        geom_ctrl_pts.push_back(p);

        // debug
//             fprintf(stderr, "t = %d geom_ctrl_pt = [%.3lf %.3lf]\n", t, geom_ctrl_pts.back().x, geom_ctrl_pts.back().y);

        vol_iter.incr_iter();
    }       // control points
    ofst.head(dom_dim) += tc.nctrl_pts;
}       // tensor products
}

// new version of geometry control points that uses anchors of tmesh
void PrepGeomCtrlPts(
    vector<vec3d>&              geom_ctrl_pts,
    Block<real_t>*              block)
{
// typing shortcuts
const auto& mfa_data        = block->mfa->geom();
const auto& tmesh           = mfa_data.tmesh;
const auto& tensor_prods    = tmesh.tensor_prods;
const auto& all_knots       = tmesh.all_knots;
const auto& all_knot_levels = tmesh.all_knot_levels;
int geom_dim                = block->mfa->geom_dim();  // number of geometry dims
int dom_dim                 = block->mfa->dom_dim;

vector<vector<float>> ctrl_pts_coords(dom_dim);
vec3d p;

// compute vectors of individual control point coordinates for the tensor product
for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
{
    const auto& tc = tensor_prods[t];

    for (auto k = 0; k < dom_dim; k++)                            // domain dimensions
    {
        for (auto j = 0; j < tc.nctrl_pts(k); j++)                      // control points
        {
            KnotIdx knot_min, idx;
            if (tc.knot_mins[k] == 0)
                knot_min = (mfa_data.p(k) + 1) / 2;
            else
                knot_min = tc.knot_mins[k];
            if (!tmesh.knot_idx_ofst(tc, knot_min, j, k, false, idx))
                throw mfa::MFAError(fmt::format("PrepGeomCtrlPts(): unable to offset knot"));

            float tsum = all_knots[k][idx];                         // odd degree, tsum is on the knot

            // odd degree, second control point from global edge is an average
            if ((mfa_data.p(k) % 2 == 1 && tc.knot_mins[k] == 0 && j == 1) ||
                    (mfa_data.p(k) % 2 == 1 && tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2))
            {
                KnotIdx idx1;
                if (tc.knot_mins[k] == 0 && j == 1)
                {
                    if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                        throw mfa::MFAError(fmt::format("PrepGeomCtrlPts(): unable to offset knot"));
                }
                else if (tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2)
                {
                    if (!tmesh.knot_idx_ofst(tc, idx, -1, k, false, idx1))
                        throw mfa::MFAError(fmt::format("PrepGeomCtrlPts(): unable to offset knot"));
                }
                tsum += all_knots[k][idx1];
                tsum /= 2.0;
            }

            if (mfa_data.p(k) % 2 == 0)                            // even degree, find center of knot span
            {
                KnotIdx idx1;
                if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                    throw mfa::MFAError(fmt::format("PrepGeomCtrlPts(): unable to offset knot"));
                tsum += all_knots[k][idx1];
                tsum /= 2.0;
            }
            ctrl_pts_coords[k].push_back(block->core_mins(k) + tsum * (block->core_maxs(k) - block->core_mins(k)));

            // debug
//                 fmt::print(stderr, "t {} k {} j {} tsum (param) {} ctrl_pts_coord {}\n",
//                         t, k, j, tsum, ctrl_pts_coords[k].back());

        }   // control points
    }   // domain dimensions
}   // tensor products

// form the tensor product of control points from the vectors of individual coordinates
VectorXi ofst = VectorXi::Zero(3);                              // offset of indices for current tensor
for (auto t = 0; t < tensor_prods.size(); t++)                  // tensor products
{
    const auto& tc   = tensor_prods[t];
    mfa::VolIterator vol_iter(tc.nctrl_pts);
    VectorXi ijk(dom_dim);
    while (!vol_iter.done())                                    // control points
    {
        vol_iter.idx_ijk(vol_iter.cur_iter(), ijk);

        if (tc.weights(vol_iter.cur_iter()) == MFA_NAW)
        {
            vol_iter.incr_iter();
            continue;
        }

        // first 3 dims stored as mesh geometry
        p.x = ctrl_pts_coords[0][ofst(0) + ijk(0)];
        if (dom_dim < 2)
            p.y = 0.0;
        else
            p.y = ctrl_pts_coords[1][ofst(1) + ijk(1)];
        if (dom_dim < 3)
            p.z = 0.0;
        else
            p.z = ctrl_pts_coords[2][ofst(2) + ijk(2)];
        geom_ctrl_pts.push_back(p);

        // debug
//             fmt::print(stderr, "t = {} geom_ctrl_pt = [{} {}]\n", t, geom_ctrl_pts.back().x, geom_ctrl_pts.back().y);

        vol_iter.incr_iter();
    }       // control points
    ofst.head(dom_dim) += tc.nctrl_pts;
}       // tensor products

}

// old version of science variable control points that averages over a window of knots
// based on Youssef's python code
// not used anymore, but kept for reference
void oldPrepSciCtrlPts(
    vector< vector <vec3d> >&   vars_ctrl_pts,
    float**&                    vars_ctrl_data,
    Block<real_t>*              block)
{
int dom_dim     = block->mfa->dom_dim;
int geom_dim    = block->mfa->geom_dim();          // number of geometry dims
int nvars           = block->mfa->nvars();                       // number of science variables
vec3d p;

vars_ctrl_pts.resize(nvars);
vars_ctrl_data = new float*[nvars];
for (size_t i = 0; i < nvars; i++)                              // science variables
{
    // typing shortcuts
    // auto& mfa_data          = block->vars[i].mfa_data;
    const auto& tmesh             = block->mfa->var(i).tmesh;
    const auto& tensor_prods      = tmesh.tensor_prods;
    const auto& all_knots         = tmesh.all_knots;
    const auto& all_knot_levels   = tmesh.all_knot_levels;

    size_t nctrl_pts = 0;
    for (auto t = 0; t < tensor_prods.size(); t++)                   // tensor products
    {
        size_t prod = 1;
        for (auto k = 0; k < dom_dim; k++)                                                        // domain dimensions
            prod *= tensor_prods[t].nctrl_pts(k);
        nctrl_pts += prod;
    }
    vars_ctrl_data[i] = new float[nctrl_pts];

    // compute vectors of individual control point coordinates for the tensor product
    vector<vector<float>> ctrl_pts_coords(geom_dim);
    for (auto t = 0; t < tensor_prods.size(); t++)              // tensor products
    {
        const TensorProduct<real_t>& tc = tensor_prods[t];
        for (auto k = 0; k < dom_dim; k++)                                                        // domain dimensions
        {
            int skip = 0;
            // starting knot in sequence for computing control point coordinate
            KnotIdx knot_min = tc.knot_mins[k];
            if (knot_min)
            {
                // skip knots at a deeper level than the tensor
                for (auto l = 0; l < block->mfa->var(i).p(k); l++)
                {
                    while (all_knot_levels[k][knot_min - l - skip] > tc.level)
                        skip++;
                }
                knot_min -= (block->mfa->var(i).p(k) - 1 + skip);
            }

            // debug
            knot_min = 0;
            fmt::print(stderr, "t {} k {} tc.knot_min {} knot_min {} skip {} tc.level {}\n",
                    t, k, tc.knot_mins[k], knot_min, skip, tc.level);

            int skip1   = skip;                                 // number of knots at a deeper level that should be skipped
            for (auto j = 0; j < tc.nctrl_pts(k); j++)              // control points
            {
                float tsum  = 0.0;
                for (auto l = 1; l < block->mfa->var(i).p(k) + 1; l++)
                {
                    // skip knots at a deeper level than the tensor
                    while (all_knot_levels[k][knot_min + j + l + skip1] > tc.level)
                        skip1++;
                    tsum += all_knots[k][knot_min + j + l + skip1];

                    // debug
                    fmt::print(stderr, "knot index {} knot value {} tsum {}\n",
                            knot_min + j + l + skip1, all_knots[k][knot_min + j + l + skip1], tsum);
                }
                tsum /= float(block->mfa->var(i).p(k));
                ctrl_pts_coords[k].push_back(block->core_mins(k) + tsum * (block->core_maxs(k) - block->core_mins(k)));

                // debug
                fmt::print(stderr, "t {} k {} j {} tsum (param) {} ctrl_pts_coord {} skip1 {}\n",
                        t, k, j, tsum, ctrl_pts_coords[k].back(), skip1);

            }   // control points
        }   // domain dimensions
    }   // tensor products

    // form the tensor product of control points from the vectors of individual coordinates
    VectorXi ofst = VectorXi::Zero(3);                              // offset of indices for current tensor
    for (auto t = 0; t < tensor_prods.size(); t++)                  // tensor products
    {
        const TensorProduct<real_t>& tc = tensor_prods[t];
        mfa::VolIterator vol_iter(tc.nctrl_pts);
        VectorXi ijk(dom_dim);
        while (!vol_iter.done())                                        // control points
        {
            vol_iter.idx_ijk(vol_iter.cur_iter(), ijk);

            if (tc.weights(vol_iter.cur_iter()) == MFA_NAW)
            {
                vol_iter.incr_iter();
                continue;
            }

            // first 3 dims stored as mesh geometry
            // control point position and optionally science variable, if the total fits in 3d
            p.x = ctrl_pts_coords[0][ofst(0) + ijk(0)];
            if (dom_dim < 2)
            {
                p.y = tc.ctrl_pts(vol_iter.cur_iter(), 0);
                p.z = 0.0;
            }
            else
            {
                p.y = ctrl_pts_coords[1][ofst(1) + ijk(1)];
                if (dom_dim < 3)
                    p.z = tc.ctrl_pts(vol_iter.cur_iter(), 0);
                else
                    p.z = ctrl_pts_coords[2][ofst(2) + ijk(2)];
            }
            vars_ctrl_pts[i].push_back(p);

            // science variable also stored as data
            // TODO this assumes a scalar variable
            vars_ctrl_data[i][vars_ctrl_pts[i].size() - 1] = tc.ctrl_pts(vol_iter.cur_iter(), 0);

            // debug
//                 fmt::print(stderr, "t {} ctrl_pt [{} {} {}]\n",
//                         t, vars_ctrl_pts[i].back().x, vars_ctrl_pts[i].back().y, vars_ctrl_data[i][vars_ctrl_pts[i].size() - 1]);

            vol_iter.incr_iter();
        }   // control points
        ofst.head(dom_dim) += tc.nctrl_pts;
    }   // tensor products
}   // science variables
}

// new version of science variable control points that uses anchors of tmesh
void PrepSciCtrlPts(
    vector<vector<vec3d>>&      vars_ctrl_pts,
    float**&                    vars_ctrl_data,
    Block<real_t>*              block)
{
int dom_dim = block->mfa->dom_dim;
int nvars   = block->mfa->nvars();
vars_ctrl_pts.resize(nvars);
vars_ctrl_data = new float*[nvars];
vec3d p;

for (size_t i = 0; i < nvars; i++)                                      // science variables
{
    // typing shortcuts
    const auto& mfa_data        = block->mfa->var(i);
    const auto& tmesh           = mfa_data.tmesh;
    const auto& tensor_prods    = tmesh.tensor_prods;
    const auto& all_knots       = tmesh.all_knots;
    const auto& all_knot_levels = tmesh.all_knot_levels;
    int var_dim                = block->mfa->var_dim(i);  // number of geometry dims

    size_t nctrl_pts = 0;
    for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
    {
        size_t prod = 1;
        for (auto k = 0; k < dom_dim; k++)
            prod *= tensor_prods[t].nctrl_pts(k);
        nctrl_pts += prod;
    }
    vars_ctrl_data[i] = new float[nctrl_pts];

    // compute vectors of individual control point coordinates for the tensor product
    vector<vector<float>> ctrl_pts_coords(dom_dim);
    for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
    {

#ifdef MFA_DEBUG_KNOT_INSERTION

            // debug: inserted control points by Boehm knot insertion
            // process these separately
            if (t == tensor_prods.size() - 1)
                continue;

#endif

        const auto& tc = tensor_prods[t];

        for (auto k = 0; k < dom_dim; k++)                            // domain dimensions
        {
            for (auto j = 0; j < tc.nctrl_pts(k); j++)                  // control points
            {
                // offset the knot to the correct control point
                KnotIdx knot_min, idx;
                if (tc.knot_mins[k] == 0)
                    knot_min = (mfa_data.p(k) + 1) / 2;
                else
                    knot_min = tc.knot_mins[k];
                if (!tmesh.knot_idx_ofst(tc, knot_min, j, k, false, idx))
                    throw mfa::MFAError(fmt::format("PrepSciCtrlPts(): unable to offset knot"));

                float tsum = all_knots[k][idx];                         // odd degree, tsum is on the knot

                // odd degree, second control point from global edge is an average
                if ((mfa_data.p(k) % 2 == 1 && tc.knot_mins[k] == 0 && j == 1) ||
                    (mfa_data.p(k) % 2 == 1 && tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2))
                {
                    KnotIdx idx1;
                    if (tc.knot_mins[k] == 0 && j == 1)
                    {
                        if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                            throw mfa::MFAError(fmt::format("PrepSciCtrlPts(): unable to offset knot"));
                    }
                    else if (tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2)
                    {
                        if (!tmesh.knot_idx_ofst(tc, idx, -1, k, false, idx1))
                            throw mfa::MFAError(fmt::format("PrepSciCtrlPts(): unable to offset knot"));
                    }
                    tsum += all_knots[k][idx1];
                    tsum /= 2.0;
                }

                if (mfa_data.p(k) % 2 == 0)                            // even degree, find center of knot span
                {
                    KnotIdx idx1;
                    if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                        throw mfa::MFAError(fmt::format("PrepSciCtrlPts(): unable to offset knot"));
                    tsum += all_knots[k][idx1];
                    tsum /= 2.0;
                }
                ctrl_pts_coords[k].push_back(block->core_mins(k) + tsum * (block->core_maxs(k) - block->core_mins(k)));

                // debug
//                     fmt::print(stderr, "t {} k {} j {} tsum (param) {} ctrl_pts_coord {}\n",
//                             t, k, j, tsum, ctrl_pts_coords[k].back());

            }   // control points
        }   // domain dimensions
    }   // tensor products

    // form the tensor product of control points from the vectors of individual coordinates
    VectorXi ofst = VectorXi::Zero(3);                                              // offset of indices for current tensor
    for (auto t = 0; t < tensor_prods.size(); t++)                            // tensor products
    {

#ifdef MFA_DEBUG_KNOT_INSERTION

            // debug: inserted control points by Boehm knot insertion
            // process these separately
            if (t == tensor_prods.size() - 1)
                continue;

#endif

        const auto& tc   = tensor_prods[t];
        mfa::VolIterator vol_iter(tc.nctrl_pts);
        VectorXi ijk(dom_dim);
        while (!vol_iter.done())                                                    // control points
        {
            vol_iter.idx_ijk(vol_iter.cur_iter(), ijk);

            if (tc.weights(vol_iter.cur_iter()) == MFA_NAW)
            {
                vol_iter.incr_iter();
                continue;
            }

            // first 3 dims stored as mesh geometry
            // control point position and optionally science variable, if the total fits in 3d
            p.x = ctrl_pts_coords[0][ofst(0) + ijk(0)];
            if (dom_dim < 2)
            {
                p.y = tc.ctrl_pts(vol_iter.cur_iter(), 0);
                p.z = 0.0;
            }
            else
            {
                p.y = ctrl_pts_coords[1][ofst(1) + ijk(1)];
                if (dom_dim < 3)
                    p.z = tc.ctrl_pts(vol_iter.cur_iter(), 0);
                else
                    p.z = ctrl_pts_coords[2][ofst(2) + ijk(2)];
            }
            vars_ctrl_pts[i].push_back(p);

            // science variable also stored as data
            vars_ctrl_data[i][vars_ctrl_pts[i].size() - 1] = tc.ctrl_pts(vol_iter.cur_iter(), 0);

            // debug
//                 fmt::print(stderr, "t {} ctrl_pt [{} {} {}]\n",
//                         t, vars_ctrl_pts[i].back().x, vars_ctrl_pts[i].back().y, vars_ctrl_data[i][vars_ctrl_pts[i].size() - 1]);

            vol_iter.incr_iter();
        }   // control points
        ofst.head(dom_dim) += tc.nctrl_pts;
    }   // tensor products
}   // science variables
}

#ifdef MFA_DEBUG_KNOT_INSERTION

// for debugging, preps inserted control points in last tensor
void PrepInsCtrlPts(
    vector< vector <vec3d> >&   vars_ctrl_pts,
    float**&                    vars_ctrl_data,
    Block<real_t>*              block)
{
int nvars = block->mfa->nvars();
vars_ctrl_pts.resize(nvars);
vars_ctrl_data = new float*[nvars];
vec3d p;

for (size_t i = 0; i < nvars; i++)                                      // science variables
{
    // typing shortcuts
    const auto& mfa_data    = block->mfa->var(i);
    const auto& tmesh             = mfa_data.tmesh;
    const auto& tensor_prods      = tmesh.tensor_prods;
    const auto& all_knots         = tmesh.all_knots;
    const auto& all_knot_levels   = tmesh.all_knot_levels;
    int geom_dim            = block->mfa->geom_dim();  // number of geometry dims
    int dom_dim             = block->mfa->dom_dim;
    const auto t                  = tensor_prods.size() - 1;
    const auto& tc                = tensor_prods[t];

    size_t nctrl_pts = 0;
    size_t prod = 1;
    for (auto k = 0; k < dom_dim; k++)
        prod *= tc.nctrl_pts(k);
    nctrl_pts += prod;
    vars_ctrl_data[i] = new float[nctrl_pts];

    // compute vectors of individual control point coordinates for the tensor product
    vector<vector<float>> ctrl_pts_coords(dom_dim);

    for (auto k = 0; k < dom_dim; k++)                            // domain dimensions
    {
        // debug: inserted control points by Boehm knot insertion
        // totally hacky for one specific case
        if (k > 0)
        {
            for (auto j = 0; j < tc.nctrl_pts(k); j++)
            {
                ctrl_pts_coords[k].push_back(block->input->domain(j * block->input->ndom_pts(k), 1));
//                     fmt::print(stderr, "ctrl_pt[{}][1] = {}\n", j, ctrl_pts_coords[k].back());
            }
            continue;
        }

        for (auto j = 0; j < tc.nctrl_pts(k); j++)                  // control points
        {
            // offset the knot to the correct control point
            KnotIdx knot_min, idx;
            if (tc.knot_mins[k] == 0)
                knot_min = (mfa_data.p(k) + 1) / 2;
            else
                knot_min = tc.knot_mins[k];
            if (!tmesh.knot_idx_ofst(tc, knot_min, j, k, false, idx))
                throw mfa::MFAError(fmt::format("PrepInsCtrlPts(): unable to offset knot"));

            float tsum = all_knots[k][idx];                         // odd degree, tsum is on the knot

            // odd degree, second control point from global edge is an average
            if ((mfa_data.p(k) % 2 == 1 && tc.knot_mins[k] == 0 && j == 1) ||
                    (mfa_data.p(k) % 2 == 1 && tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2))
            {
                KnotIdx idx1;
                if (tc.knot_mins[k] == 0 && j == 1)
                {
                    if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                        throw mfa::MFAError(fmt::format("PrepInsCtrlPts(): unable to offset knot"));
                }
                else if (tc.knot_maxs[k] == all_knots[k].size() - 1 && j == tc.nctrl_pts(k) - 2)
                {
                    if (!tmesh.knot_idx_ofst(tc, idx, -1, k, false, idx1))
                        throw mfa::MFAError(fmt::format("PrepInsCtrlPts(): unable to offset knot"));
                }
                tsum += all_knots[k][idx1];
                tsum /= 2.0;
            }

            if (mfa_data.p(k) % 2 == 0)                            // even degree, find center of knot span
            {
                KnotIdx idx1;
                if (!tmesh.knot_idx_ofst(tc, idx, 1, k, false, idx1))
                    throw mfa::MFAError(fmt::format("PrepInsCtrlPts(): unable to offset knot"));
                tsum += all_knots[k][idx1];
                tsum /= 2.0;
            }
            ctrl_pts_coords[k].push_back(block->core_mins(k) + tsum * (block->core_maxs(k) - block->core_mins(k)));

            // debug
            //                     fmt::print(stderr, "t {} k {} j {} tsum (param) {} ctrl_pts_coord {}\n",
            //                             t, k, j, tsum, ctrl_pts_coords[k].back());

        }   // control points
    }   // domain dimensions

    // form the tensor product of control points from the vectors of individual coordinates
    VectorXi ofst = VectorXi::Zero(3);                                              // offset of indices for current tensor
    mfa::VolIterator vol_iter(tc.nctrl_pts);
    VectorXi ijk(dom_dim);
    while (!vol_iter.done())                                                    // control points
    {
        vol_iter.idx_ijk(vol_iter.cur_iter(), ijk);

        if (tc.weights(vol_iter.cur_iter()) == MFA_NAW)
        {
            vol_iter.incr_iter();
            continue;
        }

        // first 3 dims stored as mesh geometry
        // control point position and optionally science variable, if the total fits in 3d
        p.x = ctrl_pts_coords[0][ofst(0) + ijk(0)];
        if (dom_dim < 2)
        {
            p.y = tc.ctrl_pts(vol_iter.cur_iter(), 0);
            p.z = 0.0;
        }
        else
        {
            p.y = ctrl_pts_coords[1][ofst(1) + ijk(1)];
            if (dom_dim < 3)
                p.z = tc.ctrl_pts(vol_iter.cur_iter(), 0);
            else
                p.z = ctrl_pts_coords[2][ofst(2) + ijk(2)];
        }
        vars_ctrl_pts[i].push_back(p);

        // science variable also stored as data
        vars_ctrl_data[i][vars_ctrl_pts[i].size() - 1] = tc.ctrl_pts(vol_iter.cur_iter(), 0);

        // debug
        //                 fmt::print(stderr, "t {} ctrl_pt [{} {} {}]\n",
        //                         t, vars_ctrl_pts[i].back().x, vars_ctrl_pts[i].back().y, vars_ctrl_data[i][vars_ctrl_pts[i].size() - 1]);

        vol_iter.incr_iter();
    }   // control points
    ofst.head(dom_dim) += tc.nctrl_pts;
}   // science variables
}

#endif

// package rendering data
void PrepRenderingData(
    vector<vec3d>&              geom_ctrl_pts,
    vector< vector <vec3d> >&   vars_ctrl_pts,
    float**&                    vars_ctrl_data,

#ifdef MFA_DEBUG_KNOT_INSERTION

    vector< vector <vec3d> >&   ins_ctrl_pts,
    float**&                    ins_ctrl_data,

#endif
    vector<vec3d>&              tensor_pts_real,
    vector<vec3d>&              tensor_pts_index,
    vector<int>&                ntensor_pts,
    Block<real_t>*              block,
    int                         sci_var)                // science variable to render geometrically for 1d and 2d domains
{
// prep control points
PrepGeomCtrlPts(geom_ctrl_pts, block);
PrepSciCtrlPts(vars_ctrl_pts, vars_ctrl_data, block);

#ifdef MFA_DEBUG_KNOT_INSERTION

PrepInsCtrlPts(ins_ctrl_pts, ins_ctrl_data, block);

#endif
// prep tmesh extents
PrepTmeshTensorExtents(tensor_pts_real, tensor_pts_index, ntensor_pts, block);
}


template<typename T>
void set_point_set(mfa::PointSet<T>*& point_set, size_t dom_dim, size_t pt_dim,
const VectorX<T>& core_min, const VectorX<T>& core_max, int upsample_factor, std::vector<T>& shrink_range_raio,Block<real_t>* block)
{
    VectorXi ndom_pts(dom_dim);
    int npts = 1;

    std::vector<T> dim_start(2*dom_dim);
    std::vector<T> modified_shrink_range_raio(2*dom_dim); //To be consisted with dim_start after rounding






    VectorXi ori_ndom_pts(dom_dim);
    auto& tc = block->mfa->var(0).tmesh.tensor_prods[0];
    VectorXi span_num = tc.nctrl_pts-block->mfa->var(0).p;
    for(int i=0;i<dom_dim;i++)
    {
        ori_ndom_pts(i) = upsample_factor * span_num(i)+1;
    }

        for(int i=0;i<2;++i)
    {
        if(shrink_range_raio[2*i+1]==1)
        {
            shrink_range_raio[2*i+1]=span_num[i];
        }
        
    }
    

    for (int i = 0; i < dom_dim; i++)
    {

        dim_start[2*i]= upsample_factor * shrink_range_raio[2*i]; //round((ori_ndom_pts[i]-1)*shrink_range_raio[2*i]);
        dim_start[2*i+1]=upsample_factor * shrink_range_raio[2*i+1];// round((ori_ndom_pts[i]-1)*shrink_range_raio[2*i+1]);

        ndom_pts(i)     =  (dim_start[2*i+1] - dim_start[2*i] + 1);
        npts    *= ndom_pts(i);

        modified_shrink_range_raio[2*i]=dim_start[2*i]/(ori_ndom_pts[i]-1);
        modified_shrink_range_raio[2*i+1]=dim_start[2*i+1]/(ori_ndom_pts[i]-1);

    }


    std::cout<<modified_shrink_range_raio[0]<<" "<<modified_shrink_range_raio[1]<<std::endl;
    std::cout<<modified_shrink_range_raio[2]<<" "<<modified_shrink_range_raio[3]<<std::endl;

    VectorXi mdims(2);
    mdims(0) = dom_dim;
    mdims(1) = pt_dim-dom_dim;
    point_set = new mfa::PointSet<T>(dom_dim, mdims, npts, ndom_pts);

    VectorX<T> d(dom_dim);               // step in domain points in each dimension
    VectorX<T> p0(dom_dim);              // starting point in each dimension
    // assign values to the domain (geometry)

    int nghost_pts;                         // number of ghost points in current dimension
    for (int i = 0; i < dom_dim; i++)
    {
        d(i) =  (core_max(i) - core_min(i)) / (ndom_pts(i) - 1) *(modified_shrink_range_raio[2*i+1]-modified_shrink_range_raio[2*i]);
        p0(i) = core_min(i)+(core_max(i) - core_min(i))*modified_shrink_range_raio[2*i];
    }

    mfa::VolIterator vol_it(ndom_pts);
    // current index of domain point in each dim, initialized to 0s
    // flattened loop over all the points in a domain
    while (!vol_it.done())
    {
        int j = (int)vol_it.cur_iter();
        // compute geometry coordinates of domain point
        for (auto i = 0; i < dom_dim; i++)
            point_set->domain(j, i) = p0(i) + vol_it.idx_dim(i) * d(i);

        vol_it.incr_iter();
    }
}



// write vtk files for initial, approximated, control points
void write_vtk_files(
        Block<real_t>* b,
        const          diy::Master::ProxyWithLink& cp,
        Block<real_t>* b2,
        const          diy::Master::ProxyWithLink& cp2,
        int            sci_var,                     // science variable to render geometrically for 1d and 2d domains
        int upsample_factor,
        std::vector<double>& shrink_range_raio,
        int ignore,
        string& output_obj_name, string& output_vtk_name)
{
    vector<vec3d>               geom_ctrl_pts;      // control points (<= 3d) in geometry
    vector < vector <vec3d> >   vars_ctrl_pts;      // control points (<= 3d) in science variables
    float**                     vars_ctrl_data;     // control point data values (4d)
#ifdef MFA_DEBUG_KNOT_INSERTION

    vector < vector <vec3d> >   ins_ctrl_pts;      // control points (<= 3d) inserted by Boehm knot insertion
    float**                     ins_ctrl_data;     // control point data values (4d)

#endif
    vector<vec3d>               tensor_pts_real;    // tmesh tensor product extents in real space
    vector<vec3d>               tensor_pts_index;   // tmesh tensor product extents in index space
    vector<int>                 ntensor_pts;        // number of tensor extent points in each dim

    // package rendering data
    PrepRenderingData(geom_ctrl_pts,
                      vars_ctrl_pts,
                      vars_ctrl_data,
#ifdef MFA_DEBUG_KNOT_INSERTION
                      ins_ctrl_pts,
                      ins_ctrl_data,
#endif
                      tensor_pts_real,
                      tensor_pts_index,
                      ntensor_pts,
                      b,
                      sci_var);

    // pad dimensions up to 3
    int dom_dim = b->mfa->dom_dim;
    int nvars   = b->mfa->nvars();
    int pt_dim = dom_dim+b->mfa->nvars();

    if(shrink_range_raio.size()<2*dom_dim)
    {
        for (auto i = (shrink_range_raio.size()>>1); i < dom_dim; i++)
        {
            shrink_range_raio.push_back(0.0);
            shrink_range_raio.push_back(1.0);
        }
    }

    // science variable settings
    int vardim          = 1;
    int centering       = 1;
    int* vardims        = new int[nvars];
    char** varnames     = new char*[nvars];
    int* centerings     = new int[nvars];
    float* vars;
    for (int i = 0; i < nvars; i++)
    {
        vardims[i]      = 1;                                // TODO; treating each variable as a scalar (for now)
        varnames[i]     = new char[256];
        centerings[i]   = 1;
        sprintf(varnames[i], "var%d", i);
    }


    if(!output_obj_name.empty())
    {
        ignore=1;
    }

    // write geometry control points
    char filename[256];


    char function_filename[256];
    
    sprintf(function_filename, output_vtk_name.c_str(), cp.gid());
    mfa::PointSet<real_t> *point_set;


            auto start_time = std::chrono::high_resolution_clock::now();
    set_point_set(point_set, dom_dim, pt_dim, b->core_mins, b->core_maxs,upsample_factor,shrink_range_raio,b);
    auto end_time = std::chrono::high_resolution_clock::now();
    // std::cout<< input->domain<<std::endl;
    std::chrono::duration<double> run_time;
    write_function_pointset_vtk(point_set,function_filename,b,b2,run_time);

            
         std::cout<<"compute value and gradient mag running time, millisecond : "<< std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time + run_time).count()/1000<<std::endl;
    

    delete[] vardims;
    for (int i = 0; i < nvars; i++)
        delete[] varnames[i];
    delete[] varnames;
    delete[] centerings;



    for (int j = 0; j < nvars; j++)
    {
        delete[] vars_ctrl_data[j];
    }
    delete[] vars_ctrl_data;

    std::cout<<"write files done"<<std::endl;
    if(output_obj_name.empty()){
        delete point_set;
    }
}



template<typename T>
void save_data(DomainArgs& d_args,string& file_name, string& input, string& output_name, int dom_dim, const VectorXi& mdims, std::vector<double>& shrink_ratio, int sci_var = -1)
{
    int total_pts = 1;
    VectorXi ndom_pts(dom_dim);
    VectorXi ori_ndom_pts(dom_dim);
    int origianl_total_points = 1;

    std::vector<int> block_dims(shrink_ratio.size());
    for (size_t i = 0; i < shrink_ratio.size(); i+=2)
    {
        block_dims[i] =round((d_args.ndom_pts[i>>1] -1) * shrink_ratio[i]);
        block_dims[i+1] = round((d_args.ndom_pts[i>>1]-1) * shrink_ratio[i+1]);

        total_pts *= block_dims[i+1]-block_dims[i] + 1;
        ndom_pts(i>>1) = block_dims[i+1]-block_dims[i] + 1;

        ori_ndom_pts(i>>1) = d_args.ndom_pts[i>>1];

        origianl_total_points*=d_args.ndom_pts[i>>1];
    }

    std::vector<T> val(origianl_total_points);
    
    // save raw data

    FILE *fd = fopen(file_name.c_str(), "r");
    assert(fd);

    std::cout<<"reading raw data from file: "<<file_name<<std::endl;

    if (!fread(&val[0], sizeof(T), origianl_total_points, fd))
    {
        perror("Error: unable to read raw file\n");
        exit(0);
    }
    

    mfa::PointSet<double>    *point_set;  
    point_set = new mfa::PointSet<double>(dom_dim, mdims, total_pts, ndom_pts);

    std::cout<<" origianl_total_points "<<origianl_total_points<<std::endl;
    std::cout<<block_dims[0]<<" "<<block_dims[1]<<" "<<block_dims[2]<<" "<<block_dims[3]<<std::endl;
    std::cout<<"ndom_pts "<<ndom_pts(0)<<" "<<ndom_pts(1)<<std::endl;
    // for (size_t i = 0; i < val.size(); i++)
    //     point_set->domain(i, 2) = val[i];

    // set domain values (just equal to i, j; ie, dx, dy = 1, 1)
    int n=0;
    if(dom_dim==2){
    for (size_t j = 0; j < (size_t)(ndom_pts(1)); j++)
        for (size_t i = 0; i < (size_t)(ndom_pts(0)); i++)
        {
            point_set->domain(n, 0) = i + block_dims[0];
            point_set->domain(n, 1) = j + block_dims[2];
            point_set->domain(n, 2) = val[ori_ndom_pts(0)*(j+block_dims[2])+i+block_dims[0]];
            n++;
        }
    }
    else if(dom_dim==3){
        for(size_t k = 0; k < (size_t)(ndom_pts(2)); k++){
            for (size_t j = 0; j < (size_t)(ndom_pts(1)); j++)
                for (size_t i = 0; i < (size_t)(ndom_pts(0)); i++)
                {
                    point_set->domain(n, 0) = i + block_dims[0];
                    point_set->domain(n, 1) = j + block_dims[2];
                    point_set->domain(n, 2) = k + block_dims[4]; //val[ori_ndom_pts(0)*(j+block_dims[2])+i+block_dims[0]];
                    n++;
                }
        }
    }

    


    char* cstr = new char[output_name.length() + 1];
    std::strcpy(cstr, output_name.c_str());
    write_pointset_vtk(point_set, cstr,sci_var);
    delete[] cstr;
}



void save_raw_data(string& file_name, string& input, string& output_name, int dom_dim, const VectorXi& mdims,std::vector<double>& shrink_ratio)
{
    std::vector<int> mdims_;
    for (int i = 0; i < mdims.size(); i++)
    {
        mdims_.push_back(mdims(i));
    }
    
    DomainArgs d_args(dom_dim, mdims_);

    if (input == "cesm")
    {
        d_args.ndom_pts.resize(2);
        d_args.ndom_pts[0]  = 3600;
        d_args.ndom_pts[1]  = 1800;

        save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }

        if (input == "vortex_street")
    {
        d_args.ndom_pts.resize(2);
        d_args.ndom_pts[0]  = 640;
        d_args.ndom_pts[1]  = 80;

        save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }

        if (input == "hurricane_isabel")
    {
        d_args.ndom_pts.resize(2);
        d_args.ndom_pts[0]  = 500;
        d_args.ndom_pts[1]  = 500;

        save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }


    if (input == "s3d")
    {
        d_args.ndom_pts.resize(3);
        d_args.ndom_pts[0]  = 704;
        d_args.ndom_pts[1]  = 540;
        d_args.ndom_pts[2]  = 550;
        save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }

    if (input == "miranda")
    {
        d_args.ndom_pts.resize(3);
        d_args.ndom_pts[0]          = 384;
        d_args.ndom_pts[1]          = 384;
        d_args.ndom_pts[2]          = 256;
        save_data<double>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }

    if (input == "nek")
    {
        d_args.ndom_pts.resize(3);
        d_args.ndom_pts[0]  = 200;
        d_args.ndom_pts[1]  = 200;
        d_args.ndom_pts[2]  = 200;
        save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);

    }

    if(input == "hurricane")
    {
        d_args.ndom_pts.resize(3);
        d_args.ndom_pts[0]          = 500;
        d_args.ndom_pts[1]          = 500;
        d_args.ndom_pts[2]     = 100;
        save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }

    if(input=="qmcpack")
    {
        d_args.ndom_pts.resize(3);
        d_args.ndom_pts[0]          = 69;
        d_args.ndom_pts[1]          = 69;
        d_args.ndom_pts[2]          = 115;
       
       save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }

     if(input=="rti")
    {
        d_args.ndom_pts.resize(3);
        d_args.ndom_pts[0]          = 144;
        d_args.ndom_pts[1]          = 256;
        d_args.ndom_pts[2]          = 256;
       
       save_data<float>(d_args,file_name, input, output_name, dom_dim, mdims,shrink_ratio);
    }
  

}



int main(int argc, char ** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);       // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;                 // equivalent of MPI_COMM_WORLD

    int                         nvars;              // number of science variables (excluding geometry)
    vector<int>                 nraw_pts;           // number of input points in each dim.
    vector<vec3d>               raw_pts;            // input raw data points (<= 3d)
    float**                     raw_data;           // input raw data values (4d)
    vector<vec3d>               geom_ctrl_pts;      // control points (<= 3d) in geometry
    vector < vector <vec3d> >   vars_ctrl_pts;      // control points (<= 3d) in science variables
    float**                     vars_ctrl_data;     // control point data values (4d)
    vector<vec3d>               approx_pts;         // aproximated data points (<= 3d)
    float**                     approx_data;        // approximated data values (4d)
    vector<vec3d>               err_pts;            // abs value error field
    string                      input  = "sinc";        // input dataset
    int                         ntest  = 0;             // number of input test points in each dim for analytical error tests
    string                      infile = "approx.mfa";  // diy input file
    string                      infile2 = "approx2.mfa";  // diy input file
    bool                        help;                   // show help
    int                         dom_dim = 2;
    int                         pt_dim = 3;        // domain and point dimensionality, respectively
    int                         sci_var = 0;            // science variable to render geometrically for 1d and 2d domains

    int upsample_factor = 1;//upsample factor for pointset based on original pointset
    // int shrink_range_raio = 1;//shrink the range of the pointset

    string input_shrink_ratio = "0-1-0-1";

    string raw_data_file;
    string output_raw_vtk = "output_raw.vtk";
    string output_vtk_name = "output.vtk";
    int ignore = 1;

    string output_obj_name;
    int output_gradient_magnitude = 0;

    // get command line arguments
    opts::Options ops;
    ops >> opts::Option('f', "infile",      infile,     " diy input file name");
    ops >> opts::Option('g', "infile",      infile2,     " diy input file name");
    ops >> opts::Option('a', "ntest",       ntest,      " number of test points in each dimension of domain (for analytical error calculation)");
    ops >> opts::Option('i', "input",       input,      " input dataset");
    ops >> opts::Option('v', "var",         sci_var,    " science variable to render geometrically for 1d and 2d domains");
    ops >> opts::Option('h', "help",        help,       " show help");
    ops >> opts::Option('u', "upsample",    upsample_factor,       " upsample factor for pointset based on original pointset");
    ops >> opts::Option('s', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-z1-z2-...\"");
    ops >> opts::Option('r', "raw data file",    raw_data_file,       " original raw data file name");
    ops >> opts::Option('o', "output raw vtk",    output_raw_vtk,       " file name of output vtk of raw file");
    ops >> opts::Option('m', "dom_dim",    dom_dim,       " domain dimensionality for raw data");
    ops >> opts::Option('d', "pt_dim",    pt_dim,       " point dimensionality for raw data");
    ops >> opts::Option('k', "ignore",    ignore,       " ignore all other files, only output MFA sampled ttk file");
    ops >> opts::Option('n', "output obj name",    output_obj_name,       " output ply file name");
    ops >> opts::Option('t', "output vtk name",    output_vtk_name,       " output vtk file name");


    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }


    std::istringstream iss(input_shrink_ratio);
    std::vector<double> shrink_ratio;
    double number;
    std::string token;
    while (std::getline(iss, token, '-')) {
        std::istringstream tokenStream(token);
        if (tokenStream >> number) {
            shrink_ratio.push_back(number);
        }
    }  


    if(2*dom_dim!=shrink_ratio.size())
    {
        for(int i=0;i<dom_dim - shrink_ratio.size()/2;i++)
        {
            shrink_ratio.push_back(0.0);
            shrink_ratio.push_back(1.0);
        }
    }

    std::cout<<"shrink ratio is: ";
    for(auto num:shrink_ratio)
    {
        std::cout<<num<<" ";
    }
    std::cout<<std::endl;

    // echo args
    fprintf(stderr, "\n--------- Input arguments ----------\n");
    cerr << "infile = " << infile << " test_points = "    << ntest <<        endl;
    if (ntest)
        cerr << "input = "          << input     << endl;
#ifdef MFA_TBB
    cerr << "threading: TBB" << endl;
#endif
#ifdef MFA_KOKKOS
    cerr << "threading: Kokkos" << endl;
#endif
#ifdef MFA_SYCL
    cerr << "threading: SYCL" << endl;
#endif
#ifdef MFA_SERIAL
    cerr << "threading: serial" << endl;
#endif
    fprintf(stderr, "-------------------------------------\n\n");



    if (!raw_data_file.empty())
    {
        VectorXi mdims(2);
        mdims(0)=dom_dim;
        mdims(1)=pt_dim-dom_dim;

        std::cout<<"save raw data"<<std::endl;
        save_raw_data(raw_data_file,input,output_raw_vtk,dom_dim,mdims,shrink_ratio);
        exit(0);
    }

    // initialize DIY
    diy::FileStorage storage("./DIY.XXXXXX");     // used for blocks to be moved out of core
    diy::Master      master(world,
            1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy);
    diy::ContiguousAssigner   assigner(world.size(), -1); // number of blocks set by read_blocks()

    std::cout<<"reading file: "<<infile<<std::endl;

    diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
    std::cout << master.size() << " blocks read from file "<< infile << "\n\n";


    diy::FileStorage storage2("./DIY.YYYYYY");     // used for blocks to be moved out of core
    diy::Master      master2(world,
            1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy);
    diy::ContiguousAssigner   assigner2(world.size(), -1); // number of blocks set by read_blocks()

    std::cout<<"reading file: "<<infile2<<std::endl;

    diy::io::read_blocks(infile2.c_str(), world, assigner2, master2, &Block<real_t>::load);
    std::cout << master2.size() << " blocks read from file "<< infile2 << "\n\n";


    // write vtk files for initial and approximated points
    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
            { 
                master2.foreach([&](Block<real_t>* b2, const diy::Master::ProxyWithLink& cp2)
                { 
                    write_vtk_files(b, cp, b2, cp2, sci_var, upsample_factor,shrink_ratio,ignore, output_obj_name,output_vtk_name); });
                    
                });
                
                
                

    // rest of the code tests analytical functions and writes those files

    if (ntest <= 0)
        exit(0);





    // master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    //         { test_and_write(b, cp, input, d_args); });
}
