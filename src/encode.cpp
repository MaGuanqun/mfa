//--------------------------------------------------------------
// encoder object
// ref: [P&T] Piegl & Tiller, The NURBS Book, 1995
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#include <mfa/mfa.hpp>
#include <mfa/encode.hpp>
#include <mfa/decode.hpp>
#include <iostream>
#include <set>

mfa::
Encoder::
Encoder(MFA& mfa_) :
    mfa(mfa_),
    p(mfa_.p),
    ndom_pts(mfa_.ndom_pts),
    nctrl_pts(mfa_.nctrl_pts),
    domain(mfa_.domain),
    params(mfa_.params),
    ctrl_pts(mfa_.ctrl_pts),
    knots(mfa_.knots),
    dom_range(mfa_.dom_range),
    po(mfa_.po),
    ko(mfa_.ko),
    knot_spans(mfa_.knot_spans)
    // DEPRECATED
//     ndone_knot_spans(mfa_.ndone_knot_spans)
{
}

// adaptive encode
void
mfa::
Encoder::
AdaptiveEncode(float err_limit)                     // maximum allowable normalized error
{
    VectorXi nnew_knots = VectorXi::Zero(p.size()); // number of new knots in each dim
    VectorXf new_knots;                             // new knots (1st dim changes fastest)
    MatrixXf approx(domain.rows(), domain.cols());  // points in approximated volume

    // fast encoding in 1-d
    for (int iter = 0; ; iter++)
    {
        fprintf(stderr, "-----\n\nEncoding iteration %d...\n", iter);
        bool done = FastEncode(nnew_knots, new_knots, err_limit);

        // debug: compute max error to see that it is decreasing
        fprintf(stderr, "\niter=%d computing max. error...\n", iter);
        mfa.Decode(approx);

        if (done)
            break;

        mfa.InsertKnots(nnew_knots, new_knots);
    }

    // TODO: slow encode full-d mfa

    // TODO: wrap everything in an outer loop until slow encode is accurate enough?
}

// fast encode using curves instead of high volume in early rounds to determine knot insertions
// returns true if done, ie, no knots are inserted
bool
mfa::
Encoder::
FastEncode(
        VectorXi& nnew_knots,                       // number of new knots in each dim
        VectorXf& new_knots,                        // new knots (1st dim changes fastest)
        float err_limit)                            // max allowable error
{
    // check and assign main quantities
    int      ndims = ndom_pts.size();               // number of domain dimensions
    VectorXi n = nctrl_pts - VectorXi::Ones(ndims); // number of control point spans in each domain dim
    VectorXi m = ndom_pts  - VectorXi::Ones(ndims); // number of input data point spans in each domain dim
    nnew_knots = VectorXi::Zero(p.size());
    new_knots.resize(0);

    // control points
    ctrl_pts.resize(mfa.tot_nctrl, domain.cols());

    for (size_t k = 0; k < ndims; k++)              // for all domain dimensions
    {
        // temporary control points for one curve
        MatrixXf temp_ctrl = MatrixXf::Zero(nctrl_pts(k), domain.cols());

        // error spans for one curve and for worst curve
        vector<int> err_spans;
        vector<int> worst_spans;
        err_spans.reserve(n(k) - p(k) + 1);
        worst_spans.reserve(n(k) - p(k) + 1);

        // maximum number of domain points with error greater than err_limit and their curves
        size_t max_nerr     =  0;
        size_t worst_curve  = -1;

        // compute the matrix N, eq. 9.66 in P&T
        // N is a matrix of (m - 1) x (n - 1) scalars that are the basis function coefficients
        //  _                                _
        // |  N_1(u[1])   ... N_{n-1}(u[1])   |
        // |     ...      ...      ...        |
        // |  N_1(u[m-1]) ... N_{n-1}(u[m-1]) |
        //  -                                -
        // TODO: N is going to be very sparse when it is large: switch to sparse representation
        // N has semibandwidth < p  nonzero entries across diagonal
        MatrixXf N = MatrixXf::Zero(m(k) - 1, n(k) - 1); // coefficients matrix

        for (int i = 1; i < m(k); i++)                  // the rows of N
        {
            int span = mfa.FindSpan(k, params(po[k] + i), ko[k]);
            assert(span - ko[k] <= n(k));            // sanity
            mfa.BasisFuns(k, params(po[k] + i), span, N, 1, n(k) - 1, i - 1, ko[k]);
        }

        // compute the product Nt x N
        // TODO: NtN is going to be very sparse when it is large: switch to sparse representation
        // NtN has semibandwidth < p + 1 nonzero entries across diagonal
        MatrixXf NtN(n(k) - 1, n(k) - 1);
        NtN = N.transpose() * N;

        // debug
        //         cerr << "k " << k << " NtN:\n" << NtN << endl;

        // R is the right hand side needed for solving NtN * P = R
        MatrixXf R(n(k) - 1, domain.cols());

        // P are the unknown interior control points and the solution to NtN * P = R
        // NtN is positive definite -> do not need pivoting
        // TODO: use a common representation for P and ctrl_pts to avoid copying
        MatrixXf P(n(k) - 1, domain.cols());

        // number of curves in this dimension
        size_t ncurves = domain.rows() / ndom_pts(k);

        size_t co = 0;                                  // starting ofst for domain curve pts in cur. dim
        size_t coo = 0;                                 // co at start of contiguous sequence
        for (size_t j = 0; j < ncurves; j++)            // for all the curves in this dimension
        {
            // compute R from input domain points
            RHS(k, N, R, ko[k], po[k], co);

            // solve for P for one curve of control points
            P = NtN.ldlt().solve(R);

            // append points from P to control points
            // TODO: any way to avoid this?
            CopyCtrl(P, n, k, co, temp_ctrl);

            // compute the error on the curve (number of input points with error > err_limit)
            int nerr = ErrorCurve(k, co, temp_ctrl, err_spans, err_limit);
            if (nerr > max_nerr)
            {
                max_nerr = nerr;
                worst_curve = j;
                worst_spans.swap(err_spans);            // shallow copy of worst_spans = err_spans
            }

            // debug
//             if (nerr)
//                 fprintf(stderr, "\ndim=%ld ncurve=%ld nerr=%d\n", k, j, nerr);

            // adjust offsets for the next curve
            if ((j + 1) % mfa.ds[k])
                co++;
            else
            {
                co = coo + mfa.ds[k] * ndom_pts(k);
                coo = co;
            }
        }                                                  // curves in this dimension

        // free R, NtN, and P
        R.resize(0, 0);
        NtN.resize(0, 0);
        P.resize(0, 0);

        // add new knots in the middle of spans with errors
        nnew_knots(k) = worst_spans.size();
        auto old_size = new_knots.size();
        new_knots.conservativeResize(old_size + worst_spans.size());    // existing values are preserved
        for (auto i = 0; i < worst_spans.size(); i++)
            new_knots(old_size + i) = (knots(ko[k] + worst_spans[i]) + knots(ko[k] + worst_spans[i] + 1)) / 2.0;

        // print progress
        fprintf(stderr, "\rdimension %ld of %d encoded", k + 1, ndims);
    }                                                      // domain dimensions

    // debug
    cerr << "\nnnew_knots:\n" << nnew_knots << endl;
    cerr << "new_knots:\n"  << new_knots  << endl;

    return(nnew_knots.sum() ? 0 : 1);
}

// approximate a NURBS hypervolume of arbitrary dimension for a given input data set
// weights are all 1 for now
// n-d version of algorithm 9.7, Piegl & Tiller (P&T) p. 422
//
// the outputs, ctrl_pts and knots, are resized by this function;  caller need not resize them
//
// There are two types of dimensionality:
// 1. The dimensionality of the NURBS tensor product (p.size())
// (1D = NURBS curve, 2D = surface, 3D = volumem 4D = hypervolume, etc.)
// 2. The dimensionality of individual domain and control points (domain.cols())
// p.size() should be < domain.cols()
void
mfa::
Encoder::
Encode()
{
    // TODO: some of these quantities mirror this in the mfa

    // check and assign main quantities
    VectorXi n;                 // number of control point spans in each domain dim
    VectorXi m;                 // number of input data point spans in each domain dim
    int      tot_nparams;       // total number of params = sum of ndom_pts over all dimensions
                                // not the total number of data points, which would be the product
    int      tot_nknots;        // total number of knots = sum of number of knots over all dims
    int      tot_nctrl;         // total number of control points
    int      ndims = ndom_pts.size();        // number of domain dimensions

    Quants(n, m, tot_nparams, tot_nknots, tot_nctrl);

    // following are counters for slicing domain and params into curves in different dimensions
    size_t cs = 1;                    // stride for domain points in curve in cur. dim
//     size_t dt = 0;                    // starting ofst for reading domain pts in curve in cur. dim.

    // control points
    ctrl_pts.resize(tot_nctrl, domain.cols());

    // 2 buffers of temporary control points
    // double buffer needed to write output curves of current dim without changing its input pts
    // temporary control points need to begin with size as many as the input domain points
    // except for the first dimension, which can be the correct number of control points
    // because the input domain points are converted to control points one dimension at a time
    // TODO: need to find a more space-efficient way
    size_t tot_ntemp_ctrl = 1;
    for (size_t k = 0; k < ndims; k++)
        tot_ntemp_ctrl *= (k == 0 ? nctrl_pts(k) : ndom_pts(k));
    MatrixXf temp_ctrl0 = MatrixXf::Zero(tot_ntemp_ctrl, domain.cols());
    MatrixXf temp_ctrl1 = MatrixXf::Zero(tot_ntemp_ctrl, domain.cols());

    VectorXi ntemp_ctrl = ndom_pts;         // current num of temp control pts in each dim

    float  max_err_val;                     // maximum solution error in final dim of all curves

    for (size_t k = 0; k < ndims; k++)      // for all domain dimensions
    {
        // TODO:
        // Investigate whether in later dimensions, when input data points are replaced by
        // control points, need new knots and params computed.
        // In the next dimension, the coordinates of the dimension didn't change,
        // but the chord length did, as control points moved away from the data points in
        // the prior dim. Need to determine how much difference it would make to recompute
        // params and knots for the new input points

        // compute the matrix N, eq. 9.66 in P&T
        // N is a matrix of (m - 1) x (n - 1) scalars that are the basis function coefficients
        //  _                                _
        // |  N_1(u[1])   ... N_{n-1}(u[1])   |
        // |     ...      ...      ...        |
        // |  N_1(u[m-1]) ... N_{n-1}(u[m-1]) |
        //  -                                -
        // TODO: N is going to be very sparse when it is large: switch to sparse representation
        // N has semibandwidth < p  nonzero entries across diagonal
        MatrixXf N = MatrixXf::Zero(m(k) - 1, n(k) - 1); // coefficients matrix

        for (int i = 1; i < m(k); i++)            // the rows of N
        {
            int span = mfa.FindSpan(k, params(po[k] + i), ko[k]);
            assert(span - ko[k] <= n(k));            // sanity
            mfa.BasisFuns(k, params(po[k] + i), span, N, 1, n(k) - 1, i - 1, ko[k]);
        }

        // debug
//         cerr << "k " << k << " N:\n" << N << endl;

        // compute the product Nt x N
        // TODO: NtN is going to be very sparse when it is large: switch to sparse representation
        // NtN has semibandwidth < p + 1 nonzero entries across diagonal
        MatrixXf NtN(n(k) - 1, n(k) - 1);
        NtN = N.transpose() * N;

        // debug
//         cerr << "k " << k << " NtN:\n" << NtN << endl;

        // R is the residual matrix needed for solving NtN * P = R
        MatrixXf R(n(k) - 1, domain.cols());

        // P are the unknown interior control points and the solution to NtN * P = R
        // NtN is positive definite -> do not need pivoting
        // TODO: use a common representation for P and ctrl_pts to avoid copying
        MatrixXf P(n(k) - 1, domain.cols());

        // number of curves in this dimension
        size_t ncurves;
        ncurves = 1;
        for (int i = 0; i < ndims; i++)
        {
            if (i < k)
                ncurves *= nctrl_pts(i);
            else if (i > k)
                ncurves *= ndom_pts(i);
            // NB: current dimension contributes no curves, hence no i == k case
        }
        // debug
        // cerr << "k: " << k << " ncurves: " << ncurves << endl;
        // cerr << "ndom_pts:\n" << ndom_pts << endl;
        // cerr << "ntemp_ctrl:\n" << ntemp_ctrl << endl;
        // if (k > 0 && k % 2 == 1) // input to odd dims is temp_ctrl0
        //     cerr << "temp_ctrl0:\n" << temp_ctrl0 << endl;
        // if (k > 0 && k % 2 == 0) // input to even dims is temp_ctrl1
        //     cerr << "temp_ctrl1:\n" << temp_ctrl1 << endl;

        size_t co = 0, to = 0;                    // starting ofst for curve & ctrl pts in cur. dim
        size_t coo = 0, too = 0;                  // co and to at start of contiguous sequence
        for (size_t j = 0; j < ncurves; j++)      // for all the curves in this dimension
        {
            // debug
            // fprintf(stderr, "j=%ld curve\n", j);

            // compute the one curve of control points
            CtrlCurve(N, NtN, R, P, n, k, co, cs, to, temp_ctrl0, temp_ctrl1);

            // adjust offsets for the next curve
            if ((j + 1) % cs)
                co++;
            else
            {
                co = coo + cs * ntemp_ctrl(k);
                coo = co;
            }
            if ((j + 1) % cs)
                to++;
            else
            {
                to = too + cs * nctrl_pts(k);
                too = to;
            }
        }                                                  // curves in this dimension

        // adjust offsets and strides for next dimension
        ntemp_ctrl(k) = nctrl_pts(k);
        cs *= ntemp_ctrl(k);

        // free R, NtN, and P
        R.resize(0, 0);
        NtN.resize(0, 0);
        P.resize(0, 0);

        // print progress
        fprintf(stderr, "\rdimension %ld of %d encoded", k + 1, ndims);

    }                                                      // domain dimensions

    fprintf(stderr,"\n");

    // debug
//     cerr << "ctrl_pts:\n" << ctrl_pts << endl;
}

// computes right hand side vector of P&T eq. 9.63 and 9.67, p. 411-412 for a curve from the
// original input domain points
void
mfa::
Encoder::
RHS(int       cur_dim,             // current dimension
    MatrixXf& N,                   // matrix of basis function coefficients
    MatrixXf& R,                   // (output) residual matrix allocated by caller
    int       ko,                  // optional index of starting knot
    int       po,                  // optional index of starting parameter
    int       co)                  // optional index of starting domain pt in current curve
{
    int n      = N.cols() + 1;               // number of control point spans
    int m      = N.rows() + 1;               // number of input data point spans

    // compute the matrix Rk for eq. 9.63 of P&T, p. 411
    MatrixXf Rk(m - 1, domain.cols());       // eigen frees MatrixX when leaving scope
    MatrixXf Nk;                             // basis coefficients for Rk[i]

    for (int k = 1; k < m; k++)
    {
        int span = mfa.FindSpan(cur_dim, params(po + k), ko);
        Nk = MatrixXf::Zero(1, n + 1);      // basis coefficients for Rk[i]
        mfa.BasisFuns(cur_dim, params(po + k), span, Nk, 0, n, 0, ko);

        Rk.row(k - 1) =
            domain.row(co + k * mfa.ds[cur_dim]) - Nk(0, 0) * domain.row(co) -
            Nk(0, n) * domain.row(co + m * mfa.ds[cur_dim]);
    }

    // compute the matrix R
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < Rk.cols(); j++)
        {
            // debug
            // fprintf(stderr, "3: i %d j %d R.rows %d Rk.rows %d\n", i, j, R.rows(), Rk.rows());
            R(i - 1, j) = (N.col(i - 1).array() * Rk.col(j).array()).sum();
        }
    }
}

// DEPRECATED; switch to the one above and confirm the answer is the same
// computes right hand side vector of P&T eq. 9.63 and 9.67, p. 411-412 for a curve from the
// original input domain points
// void
// mfa::
// Encoder::
// RHS(int       cur_dim,             // current dimension
//     MatrixXf& N,                   // matrix of basis function coefficients
//     MatrixXf& R,                   // (output) residual matrix allocated by caller
//     int       ko,                  // optional index of starting knot
//     int       po,                  // optional index of starting parameter
//     int       co,                  // optional index of starting domain pt in current curve
//     int       cs)                  // optional stride of domain pts in current curve
// {
//     int n      = N.cols() + 1;               // number of control point spans
//     int m      = N.rows() + 1;               // number of input data point spans
// 
//     // compute the matrix Rk for eq. 9.63 of P&T, p. 411
//     MatrixXf Rk(m - 1, domain.cols());       // eigen frees MatrixX when leaving scope
//     MatrixXf Nk;                             // basis coefficients for Rk[i]
// 
//     // debug
//     // cerr << "RHS domain:\n" << domain << endl;
// 
//     for (int k = 1; k < m; k++)
//     {
//         int span = mfa.FindSpan(cur_dim, params(po + k), ko);
//         Nk = MatrixXf::Zero(1, n + 1);      // basis coefficients for Rk[i]
//         mfa.BasisFuns(cur_dim, params(po + k), span, Nk, 0, n, 0, ko);
// 
//         // debug
//         // cerr << "Nk:\n" << Nk << endl;
// 
//         // debug
//         // cerr << "[" << domain.row(co + k * cs) << "] ["
//         //      << domain.row(co) << "] ["
//         //      << domain.row(co + m * cs) << "]" << endl;
// 
//         Rk.row(k - 1) =
//             domain.row(co + k * cs) - Nk(0, 0) * domain.row(co) -
//             Nk(0, n) * domain.row(co + m * cs);
//     }
// 
//     // debug
//     // cerr << "Rk:\n" << Rk << endl;
// 
//     // compute the matrix R
//     for (int i = 1; i < n; i++)
//     {
//         for (int j = 0; j < Rk.cols(); j++)
//         {
//             // debug
//             // fprintf(stderr, "3: i %d j %d R.rows %d Rk.rows %d\n", i, j, R.rows(), Rk.rows());
//             R(i - 1, j) = (N.col(i - 1).array() * Rk.col(j).array()).sum();
//         }
//     }
// }

// computes right hand side vector of P&T eq. 9.63 and 9.67, p. 411-412 for a curve from a
// new set of input points, not the default input domain
void
mfa::
Encoder::
RHS(int       cur_dim,             // current dimension
    MatrixXf& in_pts,              // input points (not the default domain stored in the mfa)
    MatrixXf& N,                   // matrix of basis function coefficients
    MatrixXf& R,                   // (output) residual matrix allocated by caller
    int       ko,                  // optional index of starting knot
    int       po,                  // optional index of starting parameter
    int       co,                  // optional index of starting input pt in current curve
    int       cs)                  // optional stride of input pts in current curve
{
    int n      = N.cols() + 1;               // number of control point spans
    int m      = N.rows() + 1;               // number of input data point spans

    // compute the matrix Rk for eq. 9.63 of P&T, p. 411
    MatrixXf Rk(m - 1, in_pts.cols());       // eigen frees MatrixX when leaving scope
    MatrixXf Nk;                             // basis coefficients for Rk[i]

    // debug
    // cerr << "RHS in_pts:\n" << in_pts << endl;

    for (int k = 1; k < m; k++)
    {
        int span = mfa.FindSpan(cur_dim, params(po + k), ko);
        Nk = MatrixXf::Zero(1, n + 1);      // basis coefficients for Rk[i]
        mfa.BasisFuns(cur_dim, params(po + k), span, Nk, 0, n, 0, ko);

        // debug
        // cerr << "Nk:\n" << Nk << endl;

        // debug
        // cerr << "[" << in_pts.row(co + k * cs) << "] ["
        //      << in_pts.row(co) << "] ["
        //      << in_pts.row(co + m * cs) << "]" << endl;

        Rk.row(k - 1) =
            in_pts.row(co + k * cs) - Nk(0, 0) * in_pts.row(co) -
            Nk(0, n) * in_pts.row(co + m * cs);
    }

    // debug
    // cerr << "Rk:\n" << Rk << endl;

    // compute the matrix R
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < Rk.cols(); j++)
        {
            // debug
            // fprintf(stderr, "3: i %d j %d R.rows %d Rk.rows %d\n", i, j, R.rows(), Rk.rows());
            R(i - 1, j) = (N.col(i - 1).array() * Rk.col(j).array()).sum();
        }
    }
}

// Checks quantities needed for approximation
void
mfa::
Encoder::
Quants(VectorXi& n,                // (output) number of control point spans in each dim
       VectorXi& m,                // (output) number of input data point spans in each dim
       int&      tot_nparams,      // (output) total number params in all dims
       int&      tot_nknots,       // (output) total number of knots in all dims
       int&      tot_nctrl)        // (output) total number of control points in all dims
{
    if (p.size() != ndom_pts.size())
    {
        fprintf(stderr, "Error: Encode() size of p must equal size of ndom_pts\n");
        exit(1);
    }
    for (size_t i = 0; i < p.size(); i++)
    {
        if (nctrl_pts(i) <= p(i))
        {
            fprintf(stderr, "Error: Encode() number of control points in dimension %ld"
                    "must be at least p + 1 for dimension %ld\n", i, i);
            exit(1);
        }
        if (nctrl_pts(i) > ndom_pts(i))
        {
            fprintf(stderr, "Error: Encode() number of control points in dimension %ld "
                    "cannot be greater than number of input data points in dimension %ld\n", i, i);
            exit(1);
        }
    }

    n.resize(p.size());
    m.resize(p.size());
    tot_nparams = 0;
    tot_nknots  = 0;
    tot_nctrl   = 1;
    for (size_t i = 0; i < p.size(); i++)
    {
        n(i)        =  nctrl_pts(i) - 1;
        m(i)        =  ndom_pts(i)  - 1;
        tot_nparams += ndom_pts(i);
        tot_nknots  += (n(i) + p(i) + 2);
        tot_nctrl   *= nctrl_pts(i);
    }
}

// append points from P to temporary control points
// init first and last control points and copy rest from solution P
// TODO: any way to avoid this copy?
// last dimension gets copied to final control points
// previous dimensions get copied to alternating double buffers
void
mfa::
Encoder::
CopyCtrl(MatrixXf& P,          // solved points for current dimension and curve
         VectorXi& n,          // number of control point spans in each dimension
         int       k,          // current dimension
         size_t    co,         // starting offset for reading domain points
         size_t    cs,         // stride for reading domain points
         size_t    to,         // starting offset for writing control points
         MatrixXf& temp_ctrl0, // first temporary control points buffer
         MatrixXf& temp_ctrl1) // second temporary control points buffer
{
    int ndims = ndom_pts.size();             // number of domain dimensions
    int nctrl_pts = n(k) + 1;                // number of control points in current dim

    // if there is only one dim, copy straight to output
    if (ndims == 1)
    {
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to, co);
        ctrl_pts.row(to) = domain.row(co);
        for (int i = 1; i < n(k); i++)
        {
            // debug
            // fprintf(stderr, "t[%ld] = p[%d]\n", to + i * cs, i - 1);
            ctrl_pts.row(to + i * cs) = P.row(i - 1);
        }
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to + n(k) * cs, co + ndom_pts(k) - 1);
        ctrl_pts.row(to + n(k) * cs) = domain.row(co + ndom_pts(k) - 1);
    }
    // first dim copied from domain to temp_ctrl0
    else if (k == 0)
    {
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to, co);
        temp_ctrl0.row(to) = domain.row(co);
        for (int i = 1; i < n(k); i++)
        {
            // debug
            // fprintf(stderr, "t[%ld] = p[%d]\n", to + i * cs, i - 1);
            temp_ctrl0.row(to + i * cs) = P.row(i - 1);
        }
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to + n(k) * cs, co + ndom_pts(k) - 1);
        temp_ctrl0.row(to + n(k) * cs) = domain.row(co + ndom_pts(k) - 1);
    }
    // even numbered dims (but not the last one) copied from temp_ctrl1 to temp_ctrl0
    else if (k % 2 == 0 && k < ndims - 1)
    {
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to, co);
        temp_ctrl0.row(to) = temp_ctrl1.row(co);
        for (int i = 1; i < n(k); i++)
        {
            // debug
            // fprintf(stderr, "t[%ld] = p[%d]\n", to + i * cs, i - 1);
            temp_ctrl0.row(to + i * cs) = P.row(i - 1);
        }
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to + n(k) * cs, co + (ndom_pts(k) - 1) * cs);
        temp_ctrl0.row(to + n(k) * cs) = temp_ctrl1.row(co + (ndom_pts(k) - 1) * cs);
    }
    // odd numbered dims (but not the last one) copied from temp_ctrl0 to temp_ctrl1
    else if (k % 2 == 1 && k < ndims - 1)
    {
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to, co);
        temp_ctrl1.row(to) = temp_ctrl0.row(co);
        for (int i = 1; i < n(k); i++)
        {
            // debug
            // fprintf(stderr, "t[%ld] = p[%d]\n", to + i * cs, i - 1);
            temp_ctrl1.row(to + i * cs) = P.row(i - 1);
        }
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to + n(k) * cs, co + (ndom_pts(k) - 1) * cs);
        temp_ctrl1.row(to + n(k) * cs) = temp_ctrl0.row(co + (ndom_pts(k) - 1) * cs);
    }
    // final dim if even is copied from temp_ctrl1 to ctrl_pts
    else if (k == ndims - 1 && k % 2 == 0)
    {
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to, co);
        ctrl_pts.row(to) = temp_ctrl1.row(co);
        for (int i = 1; i < n(k); i++)
        {
            // debug
            // fprintf(stderr, "t[%ld] = p[%d]\n", to + i * cs, i - 1);
            ctrl_pts.row(to + i * cs) = P.row(i - 1);
        }
        // debug
        // fprintf(stderr, "t[%ld] = d[%ld]\n", to + n(k) * cs, co + (ndom_pts(k) - 1) * cs);
        ctrl_pts.row(to + n(k) * cs) = temp_ctrl1.row(co + (ndom_pts(k) - 1) * cs);
    }
    // final dim if odd is copied from temp_ctrl0 to ctrl_pts
    else if (k == ndims - 1 && k % 2 == 1)
    {
        // debug
        // fprintf(stderr, "t_start[%ld] = d[%ld]\n", to, co);
        ctrl_pts.row(to) = temp_ctrl0.row(co);
        for (int i = 1; i < n(k); i++)
        {
            // debug
            // fprintf(stderr, "t[%ld] = p[%d]\n", to + i * cs, i - 1);
            ctrl_pts.row(to + i * cs) = P.row(i - 1);
        }
        // debug
        // fprintf(stderr, "t_end[%ld] = d[%ld]\n", to + n(k) * cs, co + (ndom_pts(k) - 1) * cs);
        ctrl_pts.row(to + n(k) * cs) = temp_ctrl0.row(co + (ndom_pts(k) - 1) * cs);
    }
}

// append points from P to temporary control points
// init first and last control points and copy rest from solution P
// TODO: any way to avoid this copy?
// just simple copy to one temporary buffer, no alternating double buffers
// nor copy to final control points
void
mfa::
Encoder::
CopyCtrl(MatrixXf& P,          // solved points for current dimension and curve
         VectorXi& n,          // number of control point spans in each dimension
         int       k,          // current dimension
         size_t    co,         // starting offset for reading domain points
         MatrixXf& temp_ctrl)  // temporary control points buffer
{
//     int ndims = ndom_pts.size();             // number of domain dimensions
//     int nctrl_pts = n(k) + 1;                // number of control points in current dim

    temp_ctrl.row(0) = domain.row(co);
    for (int i = 1; i < n(k); i++)
        temp_ctrl.row(i) = P.row(i - 1);
    temp_ctrl.row(n(k)) = domain.row(co + (ndom_pts(k) - 1) * mfa.ds[k]);
}

// solves for one curve of control points
void
mfa::
Encoder::
CtrlCurve(MatrixXf& N,          // basis functions for current dimension
          MatrixXf& NtN,        // N^t * N
          MatrixXf& R,          // residual matrix for current dimension and curve
          MatrixXf& P,          // solved points for current dimension and curve
          VectorXi& n,          // number of control point spans in each dimension
          size_t    k,          // current dimension
          size_t    co,         // starting ofst for reading domain pts
          size_t    cs,         // stride for reading domain points
          size_t    to,         // starting ofst for writing control pts
          MatrixXf& temp_ctrl0, // first temporary control points buffer
          MatrixXf& temp_ctrl1) // second temporary control points buffer
{
    // compute R
    // first dimension reads from domain
    // subsequent dims alternate reading temp_ctrl0 and temp_ctrl1
    // even dim reads temp_ctrl1, odd dim reads temp_ctrl0; opposite of writing order
    // because what was written in the previous dimension is read in the current one
    if (k == 0)
        RHS(k, N, R, ko[k], po[k], co);                 // input points = default domain
    else if (k % 2)
        RHS(k, temp_ctrl0, N, R, ko[k], po[k], co, cs); // input points = temp_ctrl0
    else
        RHS(k, temp_ctrl1, N, R, ko[k], po[k], co, cs); // input points = temp_ctrl1

    // solve for P
    P = NtN.ldlt().solve(R);

    // debug
    // cerr << "P:\n" << P << endl;
    // Eigen::FullPivLU<MatrixXf> lu_decomp(NtN);
    // cerr << "Rank of NtN = " << lu_decomp.rank() << endl;

    // append points from P to control points
    // TODO: any way to avoid this?
    CopyCtrl(P, n, k, co, cs, to, temp_ctrl0, temp_ctrl1);

    // debug
    // int ndims = ndom_pts.size();
    // cerr << "k " << k << " P:\n" << P << endl;
    // if (ndims == 1)
    //     cerr << "ctrl_pts:\n" << ctrl_pts << endl;
    // else if (k == 0)
    //     cerr << "temp_ctrl0:\n" << temp_ctrl0 << endl;
    // else if (k % 2 == 0 && k < ndims - 1)
    //     cerr << "temp_ctrl0:\n" << temp_ctrl0 << endl;
    // else if (k % 2 == 1 && k < ndims - 1)
    //     cerr << "temp_ctrl1:\n" << temp_ctrl1 << endl;
    // else if (k == ndims - 1)
    //     cerr << "ctrl_pts:\n" << ctrl_pts << endl;
}

// returns number of points in a curve that have error greater than err_limit
// fills err_spans with the span indices of spans that have at least one point with such error
// assumes caller reserved err_spans for max number (all) spans in the curve
// for one curve of control points
int
mfa::
Encoder::
ErrorCurve(
        size_t       k,                         // current dimension
        size_t       co,                        // starting ofst for reading domain pts
        MatrixXf&    temp_ctrl,                 // temporary control points buffer
        vector<int>& err_spans,                 // spans with error greater than err_limit
        float        err_limit)                 // max allowable error
{
    mfa::Decoder decoder(mfa);
    VectorXf cpt(domain.cols());                // decoded curve point
    int nerr = 0;                               // number of points with error greater than err_limit
    int span = p[k];                            // current knot span of the domain point being checked
    err_spans.resize(0);

    for (auto i = 0; i < ndom_pts[k]; i++)      // all domain points in the curve
    {
        while (knots(ko[k] + span + 1) < 1.0 && knots(ko[k] + span + 1) <= params(po[k] + i))
            span++;
        // debug
//         fprintf(stderr, "param=%.3f span=[%.3f %.3f]\n", params(po[k] + i), knots(ko[k] + span), knots(ko[k] + span + 1));

        decoder.CurvePt(k, params(po[k] + i), temp_ctrl, cpt);
        float err = fabs(mfa.NormalDistance(cpt, co + i * mfa.ds[k])) / dom_range;     // normalized by data range
        if (err > err_limit)
        {
            // don't duplicate spans (spans are sorted, only need to compare with back entry
            if (!err_spans.size() || span > err_spans.back())
                err_spans.push_back(span);
            nerr++;

            // debug
//             VectorXf dpt = domain.row(co + i * mfa.ds[k]);
//             cerr << "\ndomain point:\n" << dpt << endl;
//             cerr << "approx point:\n" << cpt << endl;
//             fprintf(stderr, "k=%ld i=%d co=%ld err=%.3e\n\n", k, i, co, err);
        }
    }

    return nerr;
}
