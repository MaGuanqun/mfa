//--------------------------------------------------------------
// mfa object
// ref: [P&T] Piegl & Tiller, The NURBS Book, 1995
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#include <mfa/encode.hpp>
#include <mfa/decode.hpp>

#include <Eigen/Dense>

#include <vector>
#include <iostream>

using namespace std;

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
// ------------------

mfa::
MFA::
MFA(VectorXi& p_,             // polynomial degree in each dimension
    VectorXi& ndom_pts_,      // number of input data points in each dim
    VectorXi& nctrl_pts_,     // desired number of control points in each dim
    MatrixXf& domain_,        // input data points (1st dim changes fastest)
    MatrixXf& ctrl_pts_,      // (output) control points (1st dim changes fastest)
    VectorXf& knots_) :       // (output) knots (1st dim changes fastest)
    p(p_),
    ndom_pts(ndom_pts_),
    nctrl_pts(nctrl_pts_),
    domain(domain_),
    ctrl_pts(ctrl_pts_),
    knots(knots_)
{
    // check dimensionality for sanity
    assert(p.size() < domain.cols());

    // debug
    // cerr << "domain:\n" << domain << endl;

    // total number of params = sum of ndom_pts over all dimensions
    // not the total number of data points, which would be the product
    tot_nparams = ndom_pts.sum();
    // total number of knots = sum of number of knots over all dims
    tot_nknots = 0;
    for (size_t i = 0; i < p.size(); i++)
        tot_nknots  += (nctrl_pts(i) + p(i) + 1);

    // precompute curve parameters for input points
    params.resize(tot_nparams);
    Params();

    // debug
    // cerr << "params:\n" << params << endl;

    // compute knots
    knots.resize(tot_nknots);
    Knots();

    // debug
    // cerr << "knots:\n" << knots << endl;

    // offsets and strides for knots, params, and control points in different dimensions
    // TODO: co for control points currently not used because control points are stored explicitly
    // in future, store them like params, x coords followed by y coords, ...
    // then co will be used
    ko.resize(p.size(), 0);                  // offset for knots
    po.resize(p.size(), 0);                  // offset for params
    co.resize(p.size(), 0);                  // offset for control points
    cs.resize(p.size(), 1);                  // stride for control points
    ds.resize(p.size(), 1);                  // stride for domain points
    for (size_t i = 1; i < p.size(); i++)
    {
        po[i] = po[i - 1] + ndom_pts[i - 1];
        ko[i] = ko[i - 1] + nctrl_pts[i - 1] + p[i - 1] + 1;
        co[i] = co[i - 1] * nctrl_pts[i - 1];
        ds[i] = ds[i - 1] * ndom_pts[i - 1];
    }
}

// encode
void
mfa::
MFA::
Encode()
{
    mfa::Encoder encoder(*this);
    encoder.Encode();
}

// encode
void
mfa::
MFA::
Decode(MatrixXf& approx)
{
    mfa::Decoder decoder(*this);
    decoder.Decode(approx);
}

// binary search to find the span in the knots vector containing a given parameter value
// returns span index i s.t. u is in [ knots[i], knots[i + 1] )
// NB closed interval at left and open interval at right
// except when u == knots.last(), in which case the interval is closed at both ends
// i will be in the range [p, n], where n = number of control points - 1 because there are
// p + 1 repeated knots at start and end of knot vector
// algorithm 2.1, P&T, p. 68
int
mfa::
MFA::
FindSpan(int       cur_dim,              // current dimension
         float     u,                    // parameter value
         int       ko)                   // optional starting knot to search (default = 0)
{
    if (u == knots(ko + nctrl_pts(cur_dim)))
        return ko + nctrl_pts(cur_dim) - 1;

    // binary search
    int low = p(cur_dim);
    int high = nctrl_pts(cur_dim);
    int mid = (low + high) / 2;
    while (u < knots(ko + mid) || u >= knots(ko + mid + 1))
    {
        if (u < knots(ko + mid))
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }

    return ko + mid;
}

// computes p + 1 nonvanishing basis function values [N_{span - p}, N_{span}]
// of the given parameter value
// keeps only those in the range [N_{start_n}, N_{end_n}]
// writes results in a subset of a row of N starting at index N(start_row, start_col)
// algorithm 2.2 of P&T, p. 70
// assumes N has been allocated by caller
void
mfa::
MFA::
BasisFuns(int       cur_dim,            // current dimension
          float     u,                  // parameter value
          int       span,               // index of span in the knots vector containing u
          MatrixXf& N,                  // matrix of (output) basis function values
          int       start_n,            // starting basis function N_{start_n} to compute
          int       end_n,              // ending basis function N_{end_n} to compute
          int       row,                // starting row index in N of result
          int       ko)                 // optional starting knot to search (default = 0)
{
    // init
    vector<float> scratch(p(cur_dim) + 1);            // scratchpad, same as N in P&T p. 70
    scratch[0] = 1.0;
    vector<float> left(p(cur_dim) + 1);               // temporary recurrence results
    vector<float> right(p(cur_dim) + 1);

    // fill N
    for (int j = 1; j <= p(cur_dim); j++)
    {
        left[j]  = u - knots(span + 1 - j);
        right[j] = knots(span + j) - u;
        float saved = 0.0;
        for (int r = 0; r < j; r++)
        {
            float temp = scratch[r] / (right[r + 1] + left[j - r]);
            scratch[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        scratch[j] = saved;
    }

    // debug
    // cerr << "scratch: ";
    // for (int i = 0; i < p(cur_dim) + 1; i++)
    //     cerr << scratch[i] << " ";
    // cerr << endl;

    // copy scratch to N
    for (int j = 0; j < p(cur_dim) + 1; j++)
    {
        int n_i = span - ko - p(cur_dim) + j;              // index of basis function N_{n_i}
        if (n_i >= start_n && n_i <= end_n)
        {
            int col = n_i - start_n;         // column in N where to write result
            if (col >= 0 && col < N.cols())
                N(row, col) = scratch[j];
            else
                cerr << "Note(1): BasisFuns() truncating N_" << n_i << " = " << scratch[j] <<
                    " at (" << row << ", " << col << ")" << endl;
        }
    }
}

// precompute curve parameters for input data points using the chord-length method
// n-d version of algorithm 9.3, P&T, p. 377
// params are computed along curves and averaged over all curves at same data point index i,j,k,...
// ie, resulting params for a data point i,j,k,... are same for all curves
// and params are only stored once for each dimension in row-major order (1st dim changes fastest)
// total number of params is the sum of ndom_pts over the dimensions, much less than the total
// number of data points (which would be the product)
// assumes params were allocated by caller
// TODO: investigate other schemes (domain only, normalized domain and range, etc.)
void
mfa::
MFA::
Params()
{
    float tot_dist;                    // total chord length
    VectorXf dists(ndom_pts.maxCoeff() - 1);  // chord lengths of data point spans for any dim
    params = VectorXf::Zero(params.size());
    VectorXf d;                        // current chord length

    // following are counters for slicing domain and params into curves in different dimensions
    size_t po = 0;                     // starting offset for parameters in current dim
    size_t co = 0;                     // starting offset for curves in domain in current dim
    size_t cs = 1;                     // stride for domain points in curves in current dim

    for (size_t k = 0; k < ndom_pts.size(); k++)         // for all domain dimensions
    {
        co = 0;
        size_t coo = 0;                                  // co at start of contiguous sequence
        size_t ncurves = domain.rows() / ndom_pts(k);    // number of curves in this dimension
        size_t nzero_length_curves = 0;                  // num curves with zero length
        for (size_t j = 0; j < ncurves; j++)             // for all the curves in this dimension
        {
            tot_dist = 0.0;

            // debug
            // fprintf(stderr, "1: k %d j %d po %d co %d cs %d\n", k, j, po, co, cs);

            // chord lengths
            for (size_t i = 0; i < ndom_pts(k) - 1; i++) // for all spans in this curve
            {
                // TODO: normalize domain so that dimensions they have similar scales

                // debug
                // fprintf(stderr, "  i %d co + i * cs = %d co + (i + 1) * cs = %d\n",
                //         i, co + i * cs, co + (i + 1) * cs);

                d = domain.row(co + i * cs) - domain.row(co + (i + 1) * cs);
                dists(i) = d.norm();                     // Euclidean distance (l-2 norm)
                // fprintf(stderr, "dists[%lu] = %.3f\n", i, dists[i]);
                tot_dist += dists(i);
            }

            // accumulate (sum) parameters from this curve into the params for this dim.
            if (tot_dist > 0.0)                          // skip zero length curves
            {
                params(po)                   = 0.0;      // first parameter is known
                params(po + ndom_pts(k) - 1) = 1.0;      // last parameter is known
                float prev_param             = 0.0;      // param value at previous iteration below
                for (size_t i = 0; i < ndom_pts(k) - 2; i++)
                {
                    float dfrac = dists(i) / tot_dist;
                    params(po + i + 1) += prev_param + dfrac;
                    // debug
                    // fprintf(stderr, "k %ld j %ld i %ld po %ld "
                    //         "param %.3f = prev_param %.3f + dfrac %.3f\n",
                    //         k, j, i, po, prev_param + dfrac, prev_param, dfrac);
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
            params(po + i + 1) /= (ncurves - nzero_length_curves);

        po += ndom_pts(k);
        cs *= ndom_pts(k);
    }                                                    // domain dimensions
}

// compute knots
// n-d version of eqs. 9.68, 9.69, P&T
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
// resulting knots are same for all curves and stored once for each dimension in
// row-major order (1st dim changes fastest)
// total number of knots is the sum of number of knots over the dimensions,
// much less than the product
// assumes knots were allocated by caller
void
mfa::
MFA::
Knots()
{
    // following are counters for slicing domain and params into curves in different dimensions
    size_t po = 0;                                // starting offset for params in current dim
    size_t ko = 0;                                // starting offset for knots in current dim

    for (size_t k = 0; k < p.size(); k++)         // for all domain dimensions
    {
        int nknots = nctrl_pts(k) + p(k) + 1;    // number of knots in current dim

        // in P&T, d is the ratio of number of input points (r+1) to internal knot spans (n-p+1)
        // float d = (float)(r + 1) / (n - p + 1);         // eq. 9.68, r is P&T's m
        // but I prefer d to be the ratio of input spans r to internal knot spans (n-p+1)
        float d = (float)(ndom_pts(k) - 1) / (nctrl_pts(k) - p(k));

        // compute n - p internal knots
        for (int j = 1; j <= nctrl_pts(k) - p(k) - 1; j++)    // eq. 9.69
        {
            int   i = j * d;                      // integer part of j steps of d
            float a = j * d - i;                  // fractional part of j steps of d, P&T's alpha

            // debug
            // cerr << "d " << d << " j " << j << " i " << i << " a " << a << endl;

            // when using P&T's eq. 9.68, compute knots using the following
            // knots(p + j) = (1.0 - a) * params(i - 1) + a * params(i);

            // when using my version of d, use the following
            knots(ko + p(k) + j) = (1.0 - a) * params(po + i) + a * params(po + i + 1);
        }

        // set external knots
        for (int i = 0; i < p(k) + 1; i++)
        {
            knots(ko + i) = 0.0;
            knots(ko + nknots - 1 - i) = 1.0;
        }

        po += ndom_pts(k);
        ko += nknots;
    }
}

// interpolate parameters to get parameter value for a target coordinate
//
// TODO: experiment whether this is more accurate than calling Params
// with a 1-d space of domain pts:
// min, target, and max
//
// This function currently not used, but can be useful for intersecting an MFA with an
// axis-aligned plane, which will require the parameter value corresponding to the plane coordinate
float
mfa::
MFA::
InterpolateParams(int       cur_dim,  // curent dimension
                  size_t    po,       // starting offset for params in current dim
                  size_t    ds,       // stride for domain pts in cuve in cur. dim.
                  float     coord)    // target coordinate
{
    if (coord <= domain(0, cur_dim))
        return params(po);

    if (coord >= domain((ndom_pts(cur_dim) - 1) * ds, cur_dim))
        return params(po + ndom_pts(cur_dim) - 1);

    // binary search
    int low = 0;
    int high = ndom_pts(cur_dim);
    int mid = (low + high) / 2;
    while (coord < domain((mid) * ds, cur_dim) ||
           coord >= domain((mid + 1) * ds, cur_dim))
    {
        if (coord < domain((mid) * ds, cur_dim))
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }

    // debug
    // fprintf(stderr, "binary search param po=%ld mid=%d param= %.3f\n", po, mid, params(po + mid));

    // interpolate
    // TODO: assumes the domain is ordered in increasing coordinate values
    // very dangerous!
    if (coord <= domain((mid) * ds, cur_dim) && mid > 0)
    {
        assert(coord >= domain((mid - 1) * ds, cur_dim));
        float frac = (coord - domain((mid - 1) * ds, cur_dim)) /
            (domain((mid) * ds, cur_dim) - domain((mid - 1) * ds, cur_dim));
        return params(po + mid - 1) + frac * (params(po + mid) - params(po + mid - 1));
    }
    else if (coord >= domain((mid) * ds, cur_dim) && mid < ndom_pts(cur_dim) - 1)
    {
        assert(coord <= domain((mid + 1) * ds, cur_dim));
        float frac = (coord - domain((mid) * ds, cur_dim)) /
            (domain((mid + 1) * ds, cur_dim) - domain((mid) * ds, cur_dim));
        return params(po + mid) + frac * (params(po + mid + 1) - params(po + mid));
    }
    else
        return params(po + mid, cur_dim);

    // TODO: iterate and get the param to match the target coord even closer
    // resulting coord when the param is used is within about 10^-3
}
