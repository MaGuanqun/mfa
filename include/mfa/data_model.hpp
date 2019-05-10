//--------------------------------------------------------------
// mfa data model
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------
#ifndef _DATA_MODEL_HPP
#define _DATA_MODEL_HPP

// comment out the following line for unclamped knots (single knot at each end of knot vector)
// clamped knots (repeated at ends) is the default method if no method is specified
// #define UNCLAMPED_KNOTS

// comment out the following line for domain parameterization
// domain parameterization is the default method if no method is specified
// #define CURVE_PARAMS

// comment out the following line for low-d knot insertion
// low-d is the default if no method is specified
// #define HIGH_D

// comment out the following line for applying weights to only the range dimension
// weighing the range coordinate only is the default if no method is specified
// #define WEIGH_ALL_DIMS

#include    <Eigen/Dense>
#include    <vector>
#include    <list>
#include    <iostream>

#ifndef MFA_NO_TBB
#include    <tbb/tbb.h>
using namespace tbb;
#endif

using namespace std;

using MatrixXf = Eigen::MatrixXf;
using VectorXf = Eigen::VectorXf;
using MatrixXd = Eigen::MatrixXd;
using VectorXd = Eigen::VectorXd;
using VectorXi = Eigen::VectorXi;
using ArrayXXf = Eigen::ArrayXXf;
using ArrayXXd = Eigen::ArrayXXd;
template <typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T>
using VectorX  = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using ArrayXX  = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T>
using ArrayX   = Eigen::Array<T, Eigen::Dynamic, 1>;

#include    "tmesh.hpp"

// --- data model ---
//
// using Eigen dense MartrixX to represent vectors of n-dimensional points
// rows: points; columns: point coordinates
//
// using Eigen VectorX to represent a single point
// to use a single point from a set of many points,
// explicitly copying from a row of a matrix to a vector before using the vector for math
// (and Eigen matrix row and vector are not interchangeable w/o doing an assignment,
// at least not with the current default column-major matrix ordering, not contiguous)
//
// also using Eigen VectorX to represent a set of scalars such as knots or parameters
//
// TODO: think about row/column ordering of Eigen, choose the most contiguous one
// (current default column ordering of points is not friendly to extracting a single point)
//
// TODO: think about Eigen sparse matrices
// (N and NtN matrices, used for solving for control points, are very sparse)
//
// There are two types of dimensionality:
// 1. The dimensionality of the NURBS tensor product (p.size())
// (1D = NURBS curve, 2D = surface, 3D = volumem 4D = hypervolume, etc.)
// 2. The dimensionality of individual control points (ctrl_pts.cols())
// p.size() < ctrl_pts.cols()
//
// ------------------

template <typename T>                       // float or double
struct KnotSpan
{
    VectorX<T> min_knot;                  // parameter vector of min knot in this span
    VectorX<T> max_knot;                  // parameter vector of max knot in this span
    VectorXi   min_knot_ijk;              // i,j,k indices of minimum knot in this span
    VectorXi   max_knot_ijk;              // i,j,k indices of maximum knot in this span
    VectorX<T> min_param;                 // parameter vector of minimum domain point in this span
    VectorX<T> max_param;                 // parameter vector of maximum domain point in this span
    VectorXi   min_param_ijk;             // i,j,k indices of minimum domain point in this span
    VectorXi   max_param_ijk;             // i,j,k indices of maximum domain point in this span
    int        last_split_dim;            // last dimension in which this span was subdivided
    bool       done;                      // whether the span has converged (<= error_limit everywhere)
};

namespace mfa
{
    template <typename T>                       // float or double
    class MFA_Data
    {
    public:                                     // TODO: restrict access

       VectorXi                  p;             // polynomial degree in each domain dimension
       VectorXi                  ndom_pts;      // number of input data points in each domain dim

       vector<vector<T>>         params;        // parameters for input points[dimension][index]

       Tmesh<T>                  tmesh;         // t-mesh of knots, control points, weights
       T                         range_extent;  // extent of range value of input data points
       vector< vector <size_t> > co;            // starting offset for curves in each dim

       vector<size_t>            ds;            // stride for domain points in each dim
       T                         eps;           // minimum difference considered significant
       T                         max_err;       // unnormalized absolute value of maximum error
       vector<KnotSpan <T> >     knot_spans;    // knot spans
       int                       min_dim;       // starting coordinate of this model in full-dimensional data
       int                       max_dim;       // ending coordinate of this model in full-dimensional data
       size_t                    dom_dim;       // number of domain dimensions

    public:

        // constructor for creating an mfa from input points
        MFA_Data(
                VectorXi&           p_,             // polynomial degree in each dimension
                VectorXi&           ndom_pts_,      // number of input data points in each dim
                MatrixX<T>&         domain_,        // input data points (1st dim changes fastest)
                VectorXi            nctrl_pts_,     // optional number of control points in each dim (size 0 means minimum p+1)
                int                 min_dim_,       // starting coordinate for input data
                int                 max_dim_,       // ending coordinate for input data
                T                   eps_ = 1.0e-6) :// minimum difference considered significant
            p(p_),
            ndom_pts(ndom_pts_),
            min_dim(min_dim_),
            max_dim(max_dim_),
            eps(eps_),
            tmesh(p_.size(), p_, min_dim_, max_dim_),
            dom_dim(p_.size())
        {
            // check dimensionality for sanity
            assert(dom_dim < domain_.cols());

            // max extent of input data points
            int last     = domain_.cols() - 1;
            range_extent = domain_.col(last).maxCoeff() - domain_.col(last).minCoeff();

            // set number of control points to the minimum, p + 1, if they have not been initialized
            if (!nctrl_pts_.size())
            {
                nctrl_pts_.resize(dom_dim);
                for (auto i = 0; i < dom_dim; i++)
                    nctrl_pts_(i) = p(i) + 1;
            }

            // initialize tmesh knots
            tmesh.init_knots(nctrl_pts_);

            // initialize first tensor product
            vector<size_t> knot_mins(dom_dim);
            vector<size_t> knot_maxs(dom_dim);
            for (auto i = 0; i < dom_dim; i++)
            {
                knot_mins[i] = 0;
                knot_maxs[i] = tmesh.all_knots[i].size() - 1;
            }
            tmesh.insert_tensor(knot_mins, knot_maxs);

            // stride for domain points in different dimensions
            ds.resize(dom_dim, 1);                  // stride for domain points
            for (size_t i = 1; i < dom_dim; i++)
                ds[i] = ds[i - 1] * ndom_pts[i - 1];

            // precompute curve parameters and knots for input points
            params.resize(dom_dim);

#ifdef CURVE_PARAMS
            Params(domain_);                    // params space according to the curve length (per P&T)
            Knots(tmesh);                       // knots spaced according to parameters (per P&T)
#else
            DomainParams(domain_);              // params spaced according to domain spacing
#ifndef UNCLAMPED_KNOTS
            UniformKnots(tmesh);                // knots spaced uniformly
#else
            UniformSingleKnots();                     // knots spaced uniformly with single knots at ends
#endif
#endif

            // debug
//             cerr << "Params:\n" << params << endl;
//             cerr << "Knots:\n" << knots << endl;
//             tmesh.print();

            // offsets for curve starting (domain) points in each dimension
            co.resize(dom_dim);
            for (auto k = 0; k < dom_dim; k++)
            {
                size_t ncurves  = domain_.rows() / ndom_pts(k);  // number of curves in this dimension
                size_t coo      = 0;                            // co at start of contiguous sequence
                co[k].resize(ncurves);

                co[k][0] = 0;

                for (auto j = 1; j < ncurves; j++)
                {
                    // adjust offsets for the next curve
                    if (j % ds[k])
                        co[k][j] = co[k][j - 1] + 1;
                    else
                    {
                        co[k][j] = coo + ds[k] * ndom_pts(k);
                        coo = co[k][j];
                    }
                }
            }
        }

        // constructor when reading mfa in and knowing nothing about it yet except its degree and dimensionality
        MFA_Data(
                VectorXi&           p_,             // polynomial degree in each dimension
                size_t              ntensor_prods,  // number of tensor products to allocate in tmesh
                int                 min_dim_,       // starting coordinate for input data
                int                 max_dim_,       // ending coordinate for input data
                T                   eps_ = 1.0e-6) :// minimum difference considered significant
            p(p_),
            dom_dim(p_.size()),
            min_dim(min_dim_),
            max_dim(max_dim_),
            eps(eps_),
            tmesh(p_.size(), p_, min_dim_, max_dim_, ntensor_prods) {}

        ~MFA_Data() {}

        // convert linear domain point index into (i,j,k,...) multidimensional index
        // number of dimensions is the domain dimensionality
        void idx2ijk(
                size_t    idx,                  // linear cell indx
                VectorXi& ijk)                  // i,j,k,... indices in all dimensions
        {
            if (dom_dim == 1)
            {
                ijk(0) = idx;
                return;
            }

            for (int i = 0; i < dom_dim; i++)
            {
                if (i < dom_dim - 1)
                    ijk(i) = (idx % ds[i + 1]) / ds[i];
                else
                    ijk(i) = idx / ds[i];
            }
        }

        // convert (i,j,k,...) multidimensional index into linear index into domain
        // number of dimension is the domain dimensionality
        void ijk2idx(
                VectorXi& ijk,                  // i,j,k,... indices to all dimensions
                size_t&   idx)                  // (output) linear index
        {
            idx           = 0;
            size_t stride = 1;
            for (int i = 0; i < dom_dim; i++)
            {
                idx += ijk(i) * stride;
                stride *= ndom_pts(i);
            }
        }

        // TODO: Does this need to check the level of the knots and only return the span if it's in the current tensor?
        // binary search to find the span in the knots vector containing a given parameter value
        // returns span index i s.t. u is in [ knots[i], knots[i + 1] )
        // NB closed interval at left and open interval at right
        // tmesh version
        //
        // i will be in the range [p, n], where n = number of control points - 1 because there are
        // p + 1 repeated knots at start and end of knot vector
        // algorithm 2.1, P&T, p. 68
        //
        // prints an error and returns -1 if u is not in the min,max range of knots in tensor or if levels of u and the span do not match
        int FindSpan(
                int                     cur_dim,            // current dimension
                T                       u,                  // parameter value
                const TensorProduct<T>& tensor)             // tensor product in tmesh
        {
            if (u < tmesh.all_knots[cur_dim][tensor.knot_mins[cur_dim]] ||
                    u > tmesh.all_knots[cur_dim][tensor.knot_maxs[cur_dim]])
            {
                fprintf(stderr, "FindSpan(): Asking for parameter value outside of the knot min/max of the current tensor. This should not happen.\n");
                return -1;
            }

            if (u == tmesh.all_knots[cur_dim][tensor.nctrl_pts(cur_dim)])
            {
                if (tmesh.all_knot_levels[cur_dim][tensor.nctrl_pts(cur_dim)] != tensor.level)
                {
                    fprintf(stderr, "FindSpan(): level mismatch at nctrl_pts. This should not happen.\n");
                    return -1;
                }
                return tensor.nctrl_pts(cur_dim) - 1;
            }

            // binary search
            int low = p(cur_dim);
            int high = tensor.nctrl_pts(cur_dim);
            int mid = (low + high) / 2;
            while (u < tmesh.all_knots[cur_dim][mid] || u >= tmesh.all_knots[cur_dim][mid + 1])
            {
                if (u < tmesh.all_knots[cur_dim][mid])
                    high = mid;
                else
                    low = mid;
                mid = (low + high) / 2;
            }

            while (tmesh.all_knot_levels[cur_dim][mid] > tensor.level && mid > 0)
                mid--;

            if (tmesh.all_knot_levels[cur_dim][mid] != tensor.level)
            {
                fprintf(stderr, "FindSpan(): level mismatch at mid. This should not happen.\n");
                return -1;
            }

            return mid;
        }

        // TODO: skip knots not in current level
        // computes one row of basis function values for a given parameter value
        // writes results in a row of N
        // algorithm 2.2 of P&T, p. 70
        // tmesh version
        //
        // assumes N has been allocated by caller
        void BasisFuns(
                const TensorProduct<T>& tensor,     // current tensor product
                int                     cur_dim,    // current dimension
                T                       u,          // parameter value
                int                     span,       // index of span in the knots vector containing u, relative to ko
                MatrixX<T>&             N,          // matrix of (output) basis function values
                int                     row)        // row in N of result
        {
            // init
            vector<T> scratch(p(cur_dim) + 1);                  // scratchpad, same as N in P&T p. 70
            scratch[0] = 1.0;

            // temporary recurrence results
            // left(j)  = u - knots(span + 1 - j)
            // right(j) = knots(span + j) - u
            vector<T> left(p(cur_dim) + 1);
            vector<T> right(p(cur_dim) + 1);

            // fill N
            int j_left = 1;             // j_left and j_right are like j in the loop below but skip over knots not in the right level
            int j_right = 1;
            for (int j = 1; j <= p(cur_dim); j++)
            {
                // skip knots not in current level
                while (tmesh.all_knot_levels[cur_dim][span + 1 - j_left] != tensor.level)
                {
                    j_left++;
                    assert(span + 1 - j_left >= 0);
                }
                // left[j] is u = the jth knot in the correct level to the left of span
                left[j]  = u - tmesh.all_knots[cur_dim][span + 1 - j_left];
                while (tmesh.all_knot_levels[cur_dim][span + j_right] != tensor.level)
                {
                    j_right++;
                    assert(span + j_right < tmesh.all_knot_levels[cur_dim].size());
                }
                // right[j] = the jth knot in the correct level to the right of span - u
                right[j] = tmesh.all_knots[cur_dim][span + j_right] - u;
                j_left++;
                j_right++;

                T saved = 0.0;
                for (int r = 0; r < j; r++)
                {
                    T temp = scratch[r] / (right[r + 1] + left[j - r]);
                    scratch[r] = saved + right[r + 1] * temp;
                    saved = left[j - r] * temp;
                }
                scratch[j] = saved;
            }

            // copy scratch to N
            for (int j = 0; j < p(cur_dim) + 1; j++)
                N(row, span - p(cur_dim) + j) = scratch[j];
        }

        // TODO: update to tmesh
        // computes one row of basis function values for a given parameter value
        // writes results in a row of N
        // computes first k derivatives of one row of basis function values for a given parameter value
        // output is ders, with nders + 1 rows, one for each derivative (N, N', N'', ...)
        // including origin basis functions (0-th derivatives)
        // assumes ders has been allocated by caller (nders + 1 rows, # control points cols)
        // Alg. 2.3, p. 72 of P&T
//         void DerBasisFuns(
//                 int         cur_dim,        // current dimension
//                 T           u,              // parameter value
//                 int         span,           // index of span in the knots vector containing u, relative to ko
//                 int         nders,          // number of derivatives
//                 MatrixX<T>& ders)           // output basis function derivatives
//         {
//             // matrix from p. 70 of P&T
//             // upper triangle is basis functions
//             // lower triangle is knot differences
//             MatrixX<T> ndu(p(cur_dim) + 1, p(cur_dim) + 1);
//             ndu(0, 0) = 1.0;
// 
//             // temporary recurrence results
//             // left(j)  = u - knots(span + 1 - j)
//             // right(j) = knots(span + j) - u
//             VectorX<T> left(p(cur_dim) + 1);
//             VectorX<T> right(p(cur_dim) + 1);
// 
//             // fill ndu
//             for (int j = 1; j <= p(cur_dim); j++)
//             {
//                 left(j)  = u - knots(span + ko[cur_dim] + 1 - j);
//                 right(j) = knots(span + ko[cur_dim] + j) - u;
// 
//                 T saved = 0.0;
//                 for (int r = 0; r < j; r++)
//                 {
//                     // lower triangle
//                     ndu(j, r) = right(r + 1) + left(j - r);
//                     T temp = ndu(r, j - 1) / ndu(j, r);
//                     // upper triangle
//                     ndu(r, j) = saved + right(r + 1) * temp;
//                     saved = left(j - r) * temp;
//                 }
//                 ndu(j, j) = saved;
//             }
// 
//             // two most recently computed rows a_{k,j} and a_{k-1,j}
//             MatrixX<T> a(2, p(cur_dim) + 1);
// 
//             // initialize ders and set 0-th row with the basis functions = 0-th derivatives
//             ders = MatrixX<T>::Zero(ders.rows(), ders.cols());
//             for (int j = 0; j <= p(cur_dim); j++)
//                 ders(0, span - p(cur_dim) + j) = ndu(j, p(cur_dim));
// 
//             // compute derivatives according to eq. 2.10
//             // 1st row = first derivative, 2nd row = 2nd derivative, ...
//             for (int r = 0; r <= p(cur_dim); r++)
//             {
//                 int s1, s2;                             // alternate rows in array a
//                 s1      = 0;
//                 s2      = 1;
//                 a(0, 0) = 1.0;
// 
//                 for (int k = 1; k <= nders; k++)        // over all the derivatives up to the d_th one
//                 {
//                     T d    = 0.0;
//                     int rk = r - k;
//                     int pk = p(cur_dim) - k;
// 
//                     if (r >= k)
//                     {
//                         a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
//                         d        = a(s2, 0) * ndu(rk, pk);
//                     }
// 
//                     int j1, j2;
//                     if (rk >= -1)
//                         j1 = 1;
//                     else
//                         j1 = -rk;
//                     if (r - 1 <= pk)
//                         j2 = k - 1;
//                     else
//                         j2 = p(cur_dim) - r;
// 
//                     for (int j = j1; j <= j2; j++)
//                     {
//                         a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
//                         d += a(s2, j) * ndu(rk + j, pk);
//                     }
// 
//                     if (r <= pk)
//                     {
//                         a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
//                         d += a(s2, k) * ndu(r, pk);
//                     }
// 
//                     ders(k, span - p(cur_dim) + r) = d;
//                     swap(s1, s2);
//                 }                                       // for k
//             }                                           // for r
// 
//             // multiply through by the correct factors in eq. 2.10
//             int r = p(cur_dim);
//             for (int k = 1; k <= nders; k++)
//             {
//                 ders.row(k) *= r;
//                 r *= (p(cur_dim) - k);
//             }
//         }

        // compute rational (weighted) NtN from nonrational (unweighted) N
        // ie, convert basis function coefficients to rational ones with weights
        void Rationalize(
                int         k,                      // current dimension
                VectorX<T>& weights,                // weights of control points
                MatrixX<T>& N,                      // basis function coefficients
                MatrixX<T>& NtN_rat)                // (output) rationalized Nt * N
        {
            // compute rational denominators for input points
            VectorX<T> denom(N.rows());             // rational denomoninator for param of each input point
            for (int j = 0; j < N.rows(); j++)
                denom(j) = (N.row(j).cwiseProduct(weights.transpose())).sum();

            //     cerr << "denom:\n" << denom << endl;

            // "rationalize" N and Nt
            // ie, convert their basis function coefficients to rational ones with weights
            MatrixX<T> N_rat = N;                   // need a copy because N will be reused for other curves
            for (auto i = 0; i < N.cols(); i++)
                N_rat.col(i) *= weights(i);
            for (auto j = 0; j < N.rows(); j++)
                N_rat.row(j) /= denom(j);

            // multiply rationalized Nt and N
            NtN_rat = N_rat.transpose() * N_rat;

            // debug
            //         cerr << "k " << k << " NtN:\n" << NtN << endl;
            //         cerr << " NtN_rat:\n" << NtN_rat << endl;
        }

        // signed normal distance from a point to the domain
        // uses 2-point finite differences (first order linear) method to compute gradient and normal vector
        // approximates gradient from 2 points diagonally opposite each other in all
        // domain dimensions (not from 2 points in each dimension)
        T NormalDistance(
                VectorX<T>& pt,          // point whose distance from domain is desired
                MatrixX<T>& domain,      // input data points (1st dim changes fastest)
                size_t      idx)         // index of min. corner of cell in the domain
            // that will be used to compute partial derivatives
        {
            // normal vector = [df/dx, df/dy, df/dz, ..., -1]
            // -1 is the last coordinate of the domain points, ie, the range value
            VectorX<T> normal(domain.cols());
            int      last = domain.cols() - 1;    // last coordinate of a domain pt, ie, the range value

            // convert linear idx to multidim. i,j,k... indices in each domain dimension
            VectorXi ijk(dom_dim);
            idx2ijk(idx, ijk);

            // compute i0 and i1 1d and ijk0 and ijk1 nd indices for two points in the cell in each dim.
            // even though the caller provided the minimum corner index as idx, it's
            // possible that idx is at the max end of the domain in some dimension
            // in this case we set i1 <- idx and i0 to be one less
            size_t i0, i1;                          // 1-d indices of min, max corner points
            VectorXi ijk0(dom_dim);                 // n-d ijk index of min corner
            VectorXi ijk1(dom_dim);                 // n-d ijk index of max corner
            for (int i = 0; i < dom_dim; i++)       // for all domain dimensions
            {
                // at least 2 points needed in each dimension
                // TODO: do something degenerate if not, but probably will never get to this point
                // because there will be insufficient points to encode in the first place
                assert(ndom_pts(i) >= 2);

                // two opposite corners of the cell as i,j,k coordinates
                if (ijk(i) + 1 < ndom_pts(i))
                {
                    ijk0(i) = ijk(i);
                    ijk1(i) = ijk(i) + 1;
                }
                else
                {
                    ijk0(i) = ijk(i) - 1;
                    ijk1(i) = ijk(i);
                }
            }

            // set i0 and i1 to be the 1-d indices of the corner points
            ijk2idx(ijk0, i0);
            ijk2idx(ijk1, i1);

            // compute the normal to the domain at i0 and i1
            for (int i = 0; i < dom_dim; i++)      // for all domain dimensions
                normal(i) = (domain(i1, last) - domain(i0, last)) / (domain(i1, i) - domain(i0, i));
            normal(last) = -1;
            normal /= normal.norm();

            // project distance from (pt - domain(idx)) to unit normal
            VectorX<T> dom_pt = domain.row(idx);

            // debug
            //     fprintf(stderr, "idx=%ld\n", idx);
            //     cerr << "unit normal\n" << normal << endl;
            // cerr << "point\n" << pt << endl;
            // cerr << "domain point:\n" << dom_pt << endl;
            // cerr << "pt - dom_pt:\n" << pt - dom_pt << endl;
            // fprintf(stderr, "projection = %e\n\n", normal.dot(pt - dom_pt));

            return normal.dot(pt - dom_pt);
        }

    private:

        // precompute curve parameters for input data points using the chord-length method
        // n-d version of algorithm 9.3, P&T, p. 377
        // params are computed along curves and averaged over all curves at same data point index i,j,k,...
        // ie, resulting params for a data point i,j,k,... are same for all curves
        // and params are only stored once for each dimension (1st dim params, 2nd dim params, ...)
        // total number of params is the sum of ndom_pts over the dimensions, much less than the total
        // number of data points (which would be the product)
        // assumes params were allocated by caller
        void Params(MatrixX<T>& domain)                   // input data points (1st dim changes fastest)
        {
            T          tot_dist;                          // total chord length
            VectorX<T> dists(ndom_pts.maxCoeff() - 1);    // chord lengths of data point spans for any dim
            params = VectorX<T>::Zero(params.size());
            VectorX<T> d;                                 // current chord length

            // following are counters for slicing domain and params into curves in different dimensions
            size_t co = 0;                     // starting offset for curves in domain in current dim
            size_t cs = 1;                     // stride for domain points in curves in current dim

            for (size_t k = 0; k < ndom_pts.size(); k++)         // for all domain dimensions
            {
                params[k].resize(ndom_pts(k));
                co = 0;
                size_t coo = 0;                                  // co at start of contiguous sequence
                size_t ncurves = domain.rows() / ndom_pts(k);    // number of curves in this dimension
                size_t nzero_length_curves = 0;                  // num curves with zero length
                for (size_t j = 0; j < ncurves; j++)             // for all the curves in this dimension
                {
                    tot_dist = 0.0;

                    // chord lengths
                    for (size_t i = 0; i < ndom_pts(k) - 1; i++) // for all spans in this curve
                    {
                        // TODO: normalize domain so that dimensions they have similar scales
                        d = domain.row(co + i * cs) - domain.row(co + (i + 1) * cs);
                        dists(i) = d.norm();                     // Euclidean distance (l-2 norm)
                        tot_dist += dists(i);
                    }

                    // accumulate (sum) parameters from this curve into the params for this dim.
                    if (tot_dist > 0.0)                          // skip zero length curves
                    {
                        params[k][0]                 = 0.0;      // first parameter is known
                        params[k][ndom_pts(k) - 1]   = 1.0;      // last parameter is known
                        T prev_param                 = 0.0;      // param value at previous iteration below
                        for (size_t i = 0; i < ndom_pts(k) - 2; i++)
                        {
                            T dfrac             = dists(i) / tot_dist;
                            params[k][i + 1]    += prev_param + dfrac;
                            prev_param += dfrac;
                        }
                    }
                    else
                        nzero_length_curves++;

                    if ((j + 1) % cs)
                        co++;
                    else
                    {
                        co = coo + cs * ndom_pts(k);
                        coo = co;
                    }
                }                                                // curves in this dimension

                // average the params for this dimension by dividing by the number of curves that
                // contributed to the sum (skipped zero length curves)
                for (size_t i = 0; i < ndom_pts(k) - 2; i++)
                    params[k][i + 1] /= (ncurves - nzero_length_curves);

                cs *= ndom_pts(k);
            }                                                    // domain dimensions
            // debug
            //     cerr << "params:\n" << params << endl;
        }

        // precompute parameters for input data points using domain spacing only (not length along curve)
        // params are only stored once for each dimension (1st dim params, 2nd dim params, ...)
        // total number of params is the sum of ndom_pts over the dimensions, much less than the total
        // number of data points (which would be the product)
        // assumes params were allocated by caller
        void DomainParams(MatrixX<T>& domain)                   // input data points (1st dim changes fastest)
        {
            size_t cs = 1;                                      // stride for domain points in current dim.
            for (size_t k = 0; k < dom_dim; k++)               // for all domain dimensions
            {
                params[k].resize(ndom_pts(k));
                for (size_t i = 1; i < ndom_pts(k) - 1; i++)
                    params[k][i]= fabs( (domain(cs * i, k) - domain(0, k)) /
                            (domain(cs * (ndom_pts(k) - 1), k) - domain(0, k)) );

                params[k][ndom_pts(k) - 1] = 1.0;
                cs *= ndom_pts(k);
            }                                                    // domain dimensions

            // debug
            //     cerr << "params:\n" << params << endl;
        }

        // compute knots
        // n-d version of eqs. 9.68, 9.69, P&T
        // tmesh version
        //
        // the set of knots (called U in P&T) is the set of breakpoints in the parameter space between
        // different basis functions. These are the breaks in the piecewise B-spline approximation
        //
        // nknots = n + p + 2
        // eg, for p = 3 and nctrl_pts = 7, n = nctrl_pts - 1 = 6, nknots = 11
        // let knots = {0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1}
        // there are p + 1 external knots at each end: {0, 0, 0, 0} and {1, 1, 1, 1}
        // there are n - p internal knots: {0.25, 0.5, 0.75}
        // there are n - p + 1 internal knot spans [0,0.25), [0.25, 0.5), [0.5, 0.75), [0.75, 1)
        //
        // resulting knots are same for all curves and stored once for each dimension (1st dim knots, 2nd dim, ...)
        // total number of knots is the sum of number of knots over the dimensions, much less than the product
        // assumes knots were allocated by caller
        void Knots(Tmesh<T>& tmesh)
        {
            // following are counters for slicing domain and params into curves in different dimensions
            size_t po = 0;                                          // starting offset for params in current dim
//             size_t ko = 0;                                          // starting offset for knots in current dim

            for (size_t k = 0; k < dom_dim; k++)                    // for all domain dimensions
            {
                // TODO: hard-coded for first tensor product of the tmesh
                int nctrl_pts = tmesh.tensor_prods[0].nctrl_pts(k);

                int nknots = nctrl_pts + p(k) + 1;         // number of knots in current dim

                // in P&T, d is the ratio of number of input points (r+1) to internal knot spans (n-p+1)
                //         T d = (T)(ndom_pts(k)) / (nctrl_pts - p(k));         // eq. 9.68, r is P&T's m
                // but I prefer d to be the ratio of input spans r to internal knot spans (n-p+1)
                T d = (T)(ndom_pts(k) - 1) / (nctrl_pts - p(k));

                // compute n - p internal knots
                for (int j = 1; j <= nctrl_pts - p(k) - 1; j++)  // eq. 9.69
                {
                    int   i = j * d;                                // integer part of j steps of d
                    T a = j * d - i;                                // fractional part of j steps of d, P&T's alpha

                    // when using P&T's eq. 9.68, compute knots using the following
                    //             tmesh.all_knots[k][p(k) + j] = (1.0 - a) * params(po + i - 1) + a * params(po + i);

                    // when using my version of d, use the following
                    tmesh.all_knots[k][p(k) + j] = (1.0 - a) * params(po + i) + a * params(po + i + 1);
                }

                // set external knots
                for (int i = 0; i < p(k) + 1; i++)
                {
                    tmesh.all_knots[k][i] = 0.0;
                    tmesh.all_knots[k][nknots - 1 - i] = 1.0;
                }

                po += ndom_pts(k);
            }
        }

        // compute knots
        // n-d version of uniform spacing
        // tmesh version
        //
        // nknots = n + p + 2
        // eg, for p = 3 and nctrl_pts = 7, n = nctrl_pts - 1 = 6, nknots = 11
        // let knots = {0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1}
        // there are p + 1 external knots at each end: {0, 0, 0, 0} and {1, 1, 1, 1}
        // there are n - p internal knots: {0.25, 0.5, 0.75}
        // there are n - p + 1 internal knot spans [0,0.25), [0.25, 0.5), [0.5, 0.75), [0.75, 1)
        //
        // resulting knots are same for all curves and stored once for each dimension (1st dim knots, 2nd dim, ...)
        // total number of knots is the sum of number of knots over the dimensions, much less than the product
        // assumes knots were allocated by caller
        void UniformKnots(Tmesh<T>& tmesh)
        {
            // following are counters for slicing domain and params into curves in different dimensions
            size_t po = 0;                                // starting offset for params in current dim
//             size_t ko = 0;                                // starting offset for knots in current dim

            for (size_t k = 0; k < dom_dim; k++)         // for all domain dimensions
            {
                // TODO: hard-coded for first tensor product of the tmesh
                int nctrl_pts = tmesh.tensor_prods[0].nctrl_pts(k);

                int nknots = nctrl_pts + p(k) + 1;    // number of knots in current dim

                // set p + 1 external knots at each end
                for (int i = 0; i < p(k) + 1; i++)
                {
                    tmesh.all_knots[k][i] = 0.0;
                    tmesh.all_knots[k][nknots - 1 - i] = 1.0;
                }

                // compute remaining n - p internal knots
                T step = 1.0 / (nctrl_pts - p(k));               // size of internal knot span
                for (int j = 1; j <= nctrl_pts - p(k) - 1; j++)
                    tmesh.all_knots[k][p(k) + j] = tmesh.all_knots[k][p(k) + j - 1] + step;

                po += ndom_pts(k);
            }
        }

        // TODO: change over to tmesh or deprecate, depending whether we will be splitting
        // full-dimensional knot spans and whether an index into knot spans is still needed
        // now that knots have be reorganized into vectors for each dimension and levels
        //
        // initialize knot span index
//         void KnotSpanIndex()
//         {
//             // total number of knot spans = product of number of knot spans over all dims
//             size_t int_nspans = 1;                  // number of internal (unique) spans
//             size_t all_nspans = 1;                  // total number of spans, including repeating 0s and 1s
//             for (auto i = 0; i < dom_dim; i++)
//             {
//                 int_nspans *= (nctrl_pts(i) - p(i));
//                 all_nspans *= (nctrl_pts(i) + p(i));
//             }
// 
//             knot_spans.resize(int_nspans);
// 
//             // for all knot spans, fill the KnotSpan fields
//             VectorXi ijk   = VectorXi::Zero(dom_dim);       // i,j,k of start of span
//             VectorXi p_ijk = VectorXi::Zero(dom_dim);       // i,j,k of parameter
//             size_t span_idx = 0;                            // current index into knot_spans
//             for (auto i = 0; i < all_nspans; i++)           // all knot spans (including repeated 0s and 1s)
//             {
//                 // skip repeating knot spans
//                 bool skip = false;
//                 for (auto k = 0; k < dom_dim; k++)          // dimensions
//                     if ((ijk(k) < p[k]) || ijk(k) >= nctrl_pts[k])
//                     {
//                         skip = true;
//                         break;
//                     }
// 
//                 // save knot span
//                 // TODO: may not be necessary to store all the knot span fields, but for now it is
//                 // convenient; recheck later to see which are actually used
//                 // unused ones can be computed locally below but not part of the knot span struct
//                 if (!skip)
//                 {
//                     // knot ijk
//                     knot_spans[span_idx].min_knot_ijk = ijk;
//                     knot_spans[span_idx].max_knot_ijk = ijk.array() + 1;
// 
//                     // knot values
//                     knot_spans[span_idx].min_knot.resize(dom_dim);
//                     knot_spans[span_idx].max_knot.resize(dom_dim);
//                     for (auto k = 0; k < dom_dim; k++)         // dimensions
//                     {
//                         knot_spans[span_idx].min_knot(k) = knots(ko[k] + knot_spans[span_idx].min_knot_ijk(k));
//                         knot_spans[span_idx].max_knot(k) = knots(ko[k] + knot_spans[span_idx].max_knot_ijk(k));
//                     }
// 
//                     // parameter ijk and parameter values
//                     knot_spans[span_idx].min_param.resize(dom_dim);
//                     knot_spans[span_idx].max_param.resize(dom_dim);
//                     knot_spans[span_idx].min_param_ijk.resize(dom_dim);
//                     knot_spans[span_idx].max_param_ijk.resize(dom_dim);
//                     VectorXi po_ijk = p_ijk;                    // remember starting param ijk
//                     for (auto k = 0; k < dom_dim; k++)         // dimensions in knot spans
//                     {
//                         // min param ijk and value
//                         knot_spans[span_idx].min_param_ijk(k) = p_ijk(k);
//                         knot_spans[span_idx].min_param(k)     = params(po[k] + p_ijk(k));
// 
//                         // max param ijk and value
//                         // most spans are half open [..., ...)
//                         while (params(po[k] + p_ijk(k)) < knot_spans[span_idx].max_knot(k))
//                         {
//                             knot_spans[span_idx].max_param_ijk(k) = p_ijk(k);
//                             knot_spans[span_idx].max_param(k)     = params(po[k] + p_ijk(k));
//                             p_ijk(k)++;
//                         }
//                         // the last span in each dimension is fully closed [..., ...]
//                         if (p_ijk(k) == ndom_pts(k) - 1)
//                         {
//                             knot_spans[span_idx].max_param_ijk(k) = p_ijk(k);
//                             knot_spans[span_idx].max_param(k)     = params(po[k] + p_ijk(k));
//                         }
//                     }
// 
//                     // increment param ijk
//                     for (auto k = 0; k < dom_dim; k++)     // dimension in params
//                     {
//                         if (p_ijk(k) < ndom_pts[k] - 1)
//                         {
//                             po_ijk(k) = p_ijk(k);
//                             break;
//                         }
//                         else
//                         {
//                             po_ijk(k) = 0;
//                             if (k < dom_dim - 1)
//                                 po_ijk(k + 1)++;
//                         }
//                     }
//                     p_ijk = po_ijk;
// 
//                     knot_spans[span_idx].last_split_dim = -1;
//                     knot_spans[span_idx].done           = false;
// 
//                     // debug
//                     //             cerr <<
//                     //                 "span_idx="          << span_idx                           <<
//                     //                 "\nmin_knot_ijk:\n"  << knot_spans[span_idx].min_knot_ijk  <<
//                     //                 "\nmax_knot_ijk:\n"  << knot_spans[span_idx].max_knot_ijk  <<
//                     //                 "\nmin_knot:\n"      << knot_spans[span_idx].min_knot      <<
//                     //                 "\nmax_knot:\n"      << knot_spans[span_idx].max_knot      <<
//                     //                 "\nmin_param_ijk:\n" << knot_spans[span_idx].min_param_ijk <<
//                     //                 "\nmax_param_ijk:\n" << knot_spans[span_idx].max_param_ijk <<
//                     //                 "\nmin_param:\n"     << knot_spans[span_idx].min_param     <<
//                     //                 "\nmax_param:\n"     << knot_spans[span_idx].max_param     <<
//                     //                 "\n\n"               << endl;
// 
//                     span_idx++;
//                 }                                               // !skip
// 
//                 // increment knot ijk
//                 for (auto k = 0; k < dom_dim; k++)             // dimension in knot spans
//                 {
//                     if (ijk(k) < nctrl_pts[k] + p[k] - 1)
//                     {
//                         ijk(k)++;
//                         break;
//                     }
//                     else
//                         ijk(k) = 0;
//                 }
//             }
//         }
    };
}

#endif
