// --------------------------------------------------------------
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
#include <mfa/new_knots.hpp>
#include <iostream>
#include <set>

template <typename T>                                           // float or double
mfa::
Encoder<T>::
Encoder(MFA<T>& mfa_) :
    mfa(mfa_)
{
}

// adaptive encode
template <typename T>
void
mfa::
Encoder<T>::
AdaptiveEncode(
        T   err_limit,                                          // maximum allowed normalized error
        int max_rounds)                                         // (optional) maximum number of rounds
{
    VectorXi  nnew_knots = VectorXi::Zero(mfa.p.size());        // number of new knots in each dim
    vector<T> new_knots;                                        // new knots (1st dim changes fastest)

    mfa::NewKnots<T> nk(mfa);

    // loop until no change in knots
    for (int iter = 0; ; iter++)
    {
        if (max_rounds > 0 && iter >= max_rounds)               // optional cap on number of rounds
            break;

        fprintf(stderr, "Iteration %d...\n", iter);

#ifdef HIGH_D           // high-d w/ splitting spans in the middle

        bool done = nk.NewKnots_full(nnew_knots, new_knots, err_limit, iter);

#else               // low-d w/ splitting spans in the middle

        bool done = nk.NewKnots_curve(nnew_knots, new_knots, err_limit, iter);

#endif

#if 0           // low-d w/ splitting spans at point of greatest error

        bool done = nk.NewKnots_curve1(nnew_knots, new_knots, err_limit, iter);

#endif

        // no new knots to be added
        if (done)
        {
            fprintf(stderr, "\nKnot insertion done after %d iterations; no new knots added.\n\n", iter + 1);
            break;
        }

        // check if the new knots would make the number of control points >= number of input points in any dim
        done = false;
        for (auto k = 0; k < mfa.p.size(); k++)
            if (mfa.ndom_pts(k) <= mfa.nctrl_pts(k) + nnew_knots(k))
            {
                done = true;
                break;
            }
        if (done)
        {
            fprintf(stderr, "\nKnot insertion done after %d iterations; control points would outnumber input points.\n", iter + 1);
            break;
        }
    }

    // final full encoding needed after last knot insertion above
    fprintf(stderr, "Encoding in full %ldD\n", mfa.p.size());
    Encode();
}

// linear solution of weights according to Ma and Kruth 1995 (M&K95)
// our N is M&K's B
// our Q is M&K's X_bar
// However, we skip end points that are pinned, so we solve a smaller system by 2 points in each dimension
template <typename T>
void
mfa::
Encoder<T>::
Weights(
        MatrixX<T>& Q,              // input points
        MatrixX<T>& N,              // basis function coefficients (B in M&K)
        MatrixX<T>& Nt,             // transpose of N
        MatrixX<T>& NtN,            // Nt * N
        MatrixX<T>& NtNi,           // inverse of NtN
        VectorX<T>& weights)        // output weights
{
    // TODO: decide if I need to represent all dims or only 2 dims for current curve
    // ie, we are finding weights of high-d points for a 2-d curve
    int pt_dim = mfa.domain.cols();             // dimensionality of input and control points (domain and range)
    vector<MatrixX<T>> NtQ2(pt_dim);            // temp. matrices N^T x Q^2 for each dim of points (TODO: any way to avoid?)
    vector<MatrixX<T>> NtQ2N(pt_dim);           // matrices N^T x Q^2 x N for each dim of points
    vector<MatrixX<T>> NtQ(pt_dim);             // temp. matrices N^T x Q  for each dim of points (TODO: any way to avoid?)
    vector<MatrixX<T>> NtQN(pt_dim);            // matrices N^T x Q x N for each dim of points

    // allocate matrices of NtQ, NtQ2, NtQ2N, and NtQN
    for (auto k = 0; k < pt_dim; k++)
    {
        NtQ2[k]  = MatrixX<T>::Zero(Nt.rows(),  Nt.cols());
        NtQ[k]   = MatrixX<T>::Zero(Nt.rows(),  Nt.cols());
        NtQ2N[k] = MatrixX<T>::Zero(NtN.rows(), NtN.cols());
        NtQN[k]  = MatrixX<T>::Zero(NtN.rows(), NtN.cols());
    }

    // compute NtQN and NtQ2N
    for (auto k = 0; k < pt_dim; k++)           // for all point dims
    {
        // temporary matrices NtQ and NtQ2
        for (auto i = 0; i < Nt.cols(); i++)
        {
            T dom_pt_coord = Q(i + 1, k);      // current coordinate of current input point
            NtQ[k].col(i)  = Nt.col(i) * dom_pt_coord;
            NtQ2[k].col(i) = Nt.col(i) * dom_pt_coord * dom_pt_coord;
        }
        // final matrices NtQN and NtQ2N
        NtQN[k]  = NtQ[k] * N;
        NtQ2N[k] = NtQ2[k] * N;
    }

    // compute the matrix M according to eq.3 and eq. 4 of M&K95
    MatrixX<T> M = MatrixX<T>::Zero(NtN.rows(), NtN.cols());
    for (auto k = 0; k < pt_dim; k++)           // for all point dims
        M += NtQ2N[k] - NtQN[k] * NtNi * NtQN[k];

    // compute the eigenvalues and eigenvectors of M (eq. 9 of M&K95)
    Eigen::SelfAdjointEigenSolver<MatrixX<T>> eigensolver(M);
    if (eigensolver.info() != Eigen::Success)
    {
        fprintf(stderr, "Error: Encoder::Weights(), computing eigenvalues of M failed, perhaps M is not self-adjoint?\n");
        exit(0);
    }
    // debug
//     cerr << "M:\n"            << M                          << endl;
//     cerr << "Eigenvalues:\n"  << eigensolver.eigenvalues()  << endl;
//     cerr << "Eigenvectors:\n" << eigensolver.eigenvectors() << endl;
}

#if 1

// TBB version
// ~2X faster than serial, still expensive to compute curve offsets
//
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
template <typename T>
void
mfa::
Encoder<T>::
Encode()
{
    // TODO: some of these quantities mirror this in the mfa

    // check and assign main quantities
    VectorXi n;                             // number of control point spans in each domain dim
    VectorXi m;                             // number of input data point spans in each domain dim
    int      ndims = mfa.ndom_pts.size();   // number of domain dimensions
    size_t   cs    = 1;                     // stride for input points in curve in cur. dim

    Quants(n, m);

    // control points
    mfa.ctrl_pts.resize(mfa.tot_nctrl, mfa.domain.cols());

    // 2 buffers of temporary control points
    // double buffer needed to write output curves of current dim without changing its input pts
    // temporary control points need to begin with size as many as the input domain points
    // except for the first dimension, which can be the correct number of control points
    // because the input domain points are converted to control points one dimension at a time
    // TODO: need to find a more space-efficient way
    size_t tot_ntemp_ctrl = 1;
    for (size_t k = 0; k < ndims; k++)
        tot_ntemp_ctrl *= (k == 0 ? mfa.nctrl_pts(k) : mfa.ndom_pts(k));
    MatrixX<T> temp_ctrl0 = MatrixX<T>::Zero(tot_ntemp_ctrl, mfa.domain.cols());
    MatrixX<T> temp_ctrl1 = MatrixX<T>::Zero(tot_ntemp_ctrl, mfa.domain.cols());

    VectorXi ntemp_ctrl = mfa.ndom_pts;     // current num of temp control pts in each dim

    T  max_err_val;                     // maximum solution error in final dim of all curves

    for (size_t k = 0; k < ndims; k++)      // for all domain dimensions
    {
        // number of curves in this dimension
        size_t ncurves;
        ncurves = 1;
        for (int i = 0; i < ndims; i++)
        {
            if (i < k)
                ncurves *= mfa.nctrl_pts(i);
            else if (i > k)
                ncurves *= mfa.ndom_pts(i);
            // NB: current dimension contributes no curves, hence no i == k case
        }

        // compute local version of co
        vector<size_t> co(ncurves);                     // starting curve points in current dim.
        vector<size_t> to(ncurves);                     // starting control points in current dim.
        co[0]      = 0;
        to[0]      = 0;
        size_t coo = 0;                                 // co at start of contiguous sequence
        size_t too = 0;                                 // to at start of contiguous sequence

        for (auto j = 1; j < ncurves; j++)
        {
            if (j % cs)
            {
                co[j] = co[j - 1] + 1;
                to[j] = to[j - 1] + 1;
            }
            else
            {
                co[j] = coo + cs * ntemp_ctrl(k);
                coo   = co[j];
                to[j] = too + cs * mfa.nctrl_pts(k);
                too   = to[j];
            }
            // debug
            //             fprintf(stderr, "co[%d][%d]=%ld\n", k, j, cok[j]);
        }

        // TODO:
        // Investigate whether in later dimensions, when input data points are replaced by
        // control points, need new knots and params computed.
        // In the next dimension, the coordinates of the dimension didn't change,
        // but the chord length did, as control points moved away from the data points in
        // the prior dim. Need to determine how much difference it would make to recompute
        // params and knots for the new input points
        // (moot when using domain decomposition)

        // compute the matrix N, eq. 9.66 in P&T
        // N is a matrix of (m - 1) x (n - 1) scalars that are the basis function coefficients
        //  _                                _
        // |  N_1(u[1])   ... N_{n-1}(u[1])   |
        // |     ...      ...      ...        |
        // |  N_1(u[m-1]) ... N_{n-1}(u[m-1]) |
        //  -                                -
        // TODO: N is going to be very sparse when it is large: switch to sparse representation
        // N has semibandwidth < p  nonzero entries across diagonal
        MatrixX<T> N = MatrixX<T>::Zero(m(k) - 1, n(k) - 1); // coefficients matrix

        for (int i = 1; i < m(k); i++)            // the rows of N
        {
            int span = mfa.FindSpan(k, mfa.params(mfa.po[k] + i), mfa.ko[k]) - mfa.ko[k];   // relative to ko
            assert(span <= n(k));            // sanity
            mfa.BasisFuns(k, mfa.params(mfa.po[k] + i), span, N, 1, n(k) - 1, i - 1);
        }

        // debug
//         cerr << "k " << k << " N:\n" << N << endl;

        // compute various other matrices from N
        // TODO: NtN is going to be very sparse when it is large: switch to sparse representation
        // NtN has semibandwidth < p + 1 nonzero entries across diagonal
        // TODO: I don't think I need to give the sizes, automatically allocated by the assignment
        MatrixX<T> Nt   = N.transpose();
        MatrixX<T> NtN  = Nt * N;;
        MatrixX<T> NtNi = NtN.partialPivLu().inverse();

        parallel_for (size_t(0), ncurves, [&] (size_t j)      // for all the curves in this dimension
        {
            // debug
            // fprintf(stderr, "j=%ld curve\n", j);

            // R is the right hand side needed for solving NtN * P = R
            MatrixX<T> R(n(k) - 1, mfa.domain.cols());

            // P are the unknown interior control points and the solution to NtN * P = R
            // NtN is positive definite -> do not need pivoting
            // TODO: use a common representation for P and ctrl_pts to avoid copying
            MatrixX<T> P(n(k) - 1, mfa.domain.cols());

            // compute the one curve of control points
            CtrlCurve(N, Nt, NtN, NtNi, R, P, k, co[j], cs, to[j], temp_ctrl0, temp_ctrl1);
        });                                                  // curves in this dimension

        // adjust offsets and strides for next dimension
        ntemp_ctrl(k) = mfa.nctrl_pts(k);
        cs *= ntemp_ctrl(k);

        NtN.resize(0, 0);                           // free NtN

        // print progress
        fprintf(stderr, "\rdimension %ld of %d encoded", k + 1, ndims);

    }                                                      // domain dimensions

    fprintf(stderr,"\n");

    // debug
//     cerr << "ctrl_pts:\n" << mfa.ctrl_pts << endl;
}

#else

// serial version
//
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
template <typename T>
void
mfa::
Encoder<T>::
Encode()
{
    // check and assign main quantities
    VectorXi n;                                 // number of control point spans in each domain dim
    VectorXi m;                                 // number of input data point spans in each domain dim
    int      ndims = mfa.ndom_pts.size();       // number of domain dimensions
    size_t   cs    = 1;                         // stride for domain points in curve in cur. dim

    Quants(n, m);

    // control points
    mfa.ctrl_pts.resize(mfa.tot_nctrl, mfa.domain.cols());

    // 2 buffers of temporary control points
    // double buffer needed to write output curves of current dim without changing its input pts
    // temporary control points need to begin with size as many as the input domain points
    // except for the first dimension, which can be the correct number of control points
    // because the input domain points are converted to control points one dimension at a time
    // TODO: need to find a more space-efficient way
    size_t tot_ntemp_ctrl = 1;
    for (size_t k = 0; k < ndims; k++)
        tot_ntemp_ctrl *= (k == 0 ? mfa.nctrl_pts(k) : mfa.ndom_pts(k));
    MatrixX<T> temp_ctrl0 = MatrixX<T>::Zero(tot_ntemp_ctrl, mfa.domain.cols());
    MatrixX<T> temp_ctrl1 = MatrixX<T>::Zero(tot_ntemp_ctrl, mfa.domain.cols());

    VectorXi ntemp_ctrl = mfa.ndom_pts;         // current num of temp control pts in each dim

    for (size_t k = 0; k < ndims; k++)          // for all domain dimensions
    {
        // number of curves in this dimension
        size_t ncurves;
        ncurves = 1;
        for (int i = 0; i < ndims; i++)
        {
            if (i < k)
                ncurves *= mfa.nctrl_pts(i);
            else if (i > k)
                ncurves *= mfa.ndom_pts(i);
            // NB: current dimension contributes no curves, hence no i == k case
        }

        // debug
        // cerr << "k: " << k << " ncurves: " << ncurves << endl;
        // cerr << "ndom_pts:\n" << mfa.ndom_pts << endl;
        // cerr << "ntemp_ctrl:\n" << ntemp_ctrl << endl;
        // if (k > 0 && k % 2 == 1) // input to odd dims is temp_ctrl0
        //     cerr << "temp_ctrl0:\n" << temp_ctrl0 << endl;
        // if (k > 0 && k % 2 == 0) // input to even dims is temp_ctrl1
        //     cerr << "temp_ctrl1:\n" << temp_ctrl1 << endl;

        // compute local version of co
        vector<size_t> co(ncurves);                     // starting curve points in current dim.
        vector<size_t> to(ncurves);                     // starting control points in current dim.
        co[0]      = 0;
        to[0]      = 0;
        size_t coo = 0;                                 // co at start of contiguous sequence
        size_t too = 0;                                 // to at start of contiguous sequence

        for (auto j = 1; j < ncurves; j++)
        {
            if (j % cs)
            {
                co[j] = co[j - 1] + 1;
                to[j] = to[j - 1] + 1;
            }
            else
            {
                co[j] = coo + cs * ntemp_ctrl(k);
                coo   = co[j];
                to[j] = too + cs * mfa.nctrl_pts(k);
                too   = to[j];
            }
            // debug
            //             fprintf(stderr, "co[%d][%d]=%ld\n", k, j, cok[j]);
        }

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
        MatrixX<T> N = MatrixX<T>::Zero(m(k) - 1, n(k) - 1); // coefficients matrix

        for (int i = 1; i < m(k); i++)            // the rows of N
        {
            int span = mfa.FindSpan(k, mfa.params(mfa.po[k] + i), mfa.ko[k]) - mfa.ko[k];  // relative to ko
            assert(span <= n(k));            // sanity
            mfa.BasisFuns(k, mfa.params(mfa.po[k] + i), span, N, 1, n(k) - 1, i - 1);
        }

        // debug
//         cerr << "k " << k << " N:\n" << N << endl;
//         for (auto i = 0; i < N.rows(); i++)
//             cerr << "row " << i << " sum (should be 1.0): " << N.row(i).sum() << endl;

        // compute various other matrices from N
        // TODO: NtN is going to be very sparse when it is large: switch to sparse representation
        // NtN has semibandwidth < p + 1 nonzero entries across diagonal
        // TODO: I don't think I need to give the sizes, automatically allocated by the assignment
        MatrixX<T> Nt   = N.transpose();
        MatrixX<T> NtN  = Nt * N;;
        MatrixX<T> NtNi = NtN.partialPivLu().inverse();

        // R is the residual matrix needed for solving NtN * P = R
        MatrixX<T> R(n(k) - 1, mfa.domain.cols());

        // P are the unknown interior control points and the solution to NtN * P = R
        // NtN is positive definite -> do not need pivoting
        // TODO: use a common representation for P and ctrl_pts to avoid copying
        MatrixX<T> P(n(k) - 1, mfa.domain.cols());

        // encode curves in this dimension
        for (size_t j = 0; j < ncurves; j++)
        {
            // print progress
            if (j > 0 && j > 100 && j % (ncurves / 100) == 0)
                fprintf(stderr, "\r dimension %ld: %.0f %% encoded (%ld out of %ld curves)",
                        k, (T)j / (T)ncurves * 100, j, ncurves);

            // compute the one curve of control points
            CtrlCurve(N, Nt, NtN, NtNi, R, P, k, co[j], cs, to[j], temp_ctrl0, temp_ctrl1);
        }

        // adjust offsets and strides for next dimension
        ntemp_ctrl(k) = mfa.nctrl_pts(k);
        cs *= ntemp_ctrl(k);

        // free R, NtN, and P
        R.resize(0, 0);
        NtN.resize(0, 0);
        P.resize(0, 0);

        // print progress
        fprintf(stderr, "\33[2K\rdimension %ld of %d encoded\n", k + 1, ndims);

    }                                                      // domain dimensions

    fprintf(stderr,"\n");

    // debug
//     cerr << "ctrl_pts:\n" << mfa.ctrl_pts << endl;
}

#endif

// computes right hand side vector of P&T eq. 9.63 and 9.67, p. 411-412 for a curve from the
// original input domain points
//
// includes multiplication by weights, so that R = (N^TQ)Nw
// where:
// N is the matrix of basis function coefficients
// Q is a diagonal matrix of input points on the diagonals and 0 elsewhere (excluding first and last points)
// w is the weights vector, excluding first and last weights
// R is column vector of n - 1 elements, each element multiple coordinates of the input points
// ie, a matrix of n - 1 rows and control point dims columns
template <typename T>
void
mfa::
Encoder<T>::
RHS(
        int         cur_dim,             // current dimension
        MatrixX<T>& N,                   // matrix of basis function coefficients
        MatrixX<T>& R,                   // (output) residual matrix allocated by caller
        VectorX<T>& weights,             // precomputed weights for n + 1 control points on this curve
        int         ko,                  // index of starting knot
        int         po,                  // index of starting parameter
        int         co)                  // index of starting domain pt in current curve
{
    int n      = N.cols() + 1;               // number of control point spans
    int m      = N.rows() + 1;               // number of input data point spans

    // compute the matrix Rk for eq. 9.63 of P&T, p. 411
    MatrixX<T> Rk(m - 1, mfa.domain.cols());   // eigen frees MatrixX when leaving scope
    MatrixX<T> Nk;                             // basis coefficients for Rk[i]
    VectorX<T> denom(m - 1);                   // rational denomoninator for param of each input point excl. ends

    for (int k = 1; k < m; k++)
    {
        int span = mfa.FindSpan(cur_dim, mfa.params(po + k), ko) - ko;     // relative to ko
        Nk = MatrixX<T>::Zero(1, n + 1);      // basis coefficients for Rk[i]
        mfa.BasisFuns(cur_dim, mfa.params(po + k), span, Nk, 0, n, 0);

        denom(k - 1) = (Nk.row(0).cwiseProduct(weights.transpose())).sum();

        Rk.row(k - 1) =
            mfa.domain.row(co + k * mfa.ds[cur_dim]) -
            Nk(0, 0) * weights(0) / denom(k - 1) * mfa.domain.row(co) -
            Nk(0, n) * weights(0) / denom(k - 1) * mfa.domain.row(co + m * mfa.ds[cur_dim]);
    }

    // compute the matrix R
    for (int i = 1; i < n; i++)
        for (int j = 0; j < Rk.cols(); j++)
            // using array() for element-wise multiplication, which is what we want (not matrix mult.)
            R(i - 1, j) =
                (N.col(i - 1).array() *             // ith basis functions for input pts
                 weights(i) / denom.array() *       // rationalized
                 Rk.col(j).array()).sum();          // input points

    // debug
//     cerr << "R:\n" << R << endl;
}

// computes right hand side vector of P&T eq. 9.63 and 9.67, p. 411-412 for a curve from a
// new set of input points, not the default input domain
//
// includes multiplication by weights, so that R = (N^TQ)Nw
// where:
// N is the matrix of basis function coefficients
// Q is a diagonal matrix of input points on the diagonals and 0 elsewhere (excluding first and last points)
// w is the weights vector, excluding first and last weights
// R is column vector of n - 1 elements, each element multiple coordinates of the input points
// ie, a matrix of n - 1 rows and control point dims columns
template <typename T>
void
mfa::
Encoder<T>::
RHS(
        int         cur_dim,             // current dimension
        MatrixX<T>& in_pts,              // input points (not the default domain stored in the mfa)
        MatrixX<T>& N,                   // matrix of basis function coefficients
        MatrixX<T>& R,                   // (output) residual matrix allocated by caller
        VectorX<T>& weights,             // precomputed weights for n + 1 control points on this curve
        int         ko,                  // index of starting knot
        int         po,                  // index of starting parameter
        int         co,                  // index of starting input pt in current curve
        int         cs)                  // stride of input pts in current curve
{
    int n      = N.cols() + 1;               // number of control point spans
    int m      = N.rows() + 1;               // number of input data point spans

    // compute the matrix Rk for eq. 9.63 of P&T, p. 411
    MatrixX<T> Rk(m - 1, in_pts.cols());       // eigen frees MatrixX when leaving scope
    MatrixX<T> Nk;                             // basis coefficients for Rk[i]
    VectorX<T> denom(m - 1);                   // rational denomoninator for param of each input point excl. ends

    for (int k = 1; k < m; k++)
    {
        int span = mfa.FindSpan(cur_dim, mfa.params(po + k), ko) - ko;
        Nk = MatrixX<T>::Zero(1, n + 1);      // basis coefficients for Rk[i]
        mfa.BasisFuns(cur_dim, mfa.params(po + k), span, Nk, 0, n, 0);

        denom(k - 1) = (Nk.row(0).cwiseProduct(weights.transpose())).sum();

        Rk.row(k - 1) =
            in_pts.row(co + k * cs) -
            Nk(0, 0) * weights(0) / denom(k - 1) * in_pts.row(co) -
            Nk(0, n) * weights(0) / denom(k - 1) * in_pts.row(co + m * cs);
    }

    // compute the matrix R
    for (int i = 1; i < n; i++)
        for (int j = 0; j < Rk.cols(); j++)
            R(i - 1, j) =
                (N.col(i - 1).array() *             // ith basis functions for input pts
                 weights(i) / denom.array() *       // rationalized
                 Rk.col(j).array()).sum();          // input points

    // debug
//     cerr << "R:\n" << R << endl;
}

// computes weighted right hand side vector of M&K'95 eq. 5 premultiplied by N^T, ie R = N^T * Q * N * w
// input is a curve from the original input domain points
//
// where:
// N is the matrix of basis function coefficients
// Q is a diagonal matrix of input points on the diagonals and 0 elsewhere (excluding first and last points)
// w is the weights vector, excluding first and last weights
// R is column vector of n - 1 elements, each element multiple coordinates of the input points
// ie, a matrix of n - 1 rows and control point dims columns
template <typename T>
void
mfa::
Encoder<T>::
RHS(
        int         cur_dim,             // current dimension
        MatrixX<T>& N,                   // matrix of basis function coefficients
        MatrixX<T>& Nt,                  // transpose of N
        MatrixX<T>& R,                   // (output) residual matrix allocated by caller
        VectorX<T>& weights,             // precomputed weights for n + 1 control points on this curve
        int         ko,                  // index of starting knot
        int         po,                  // index of starting parameter
        int         co)                  // index of starting domain pt in current curve
{
    // TODO: decide if I need to represent all dims or only 2 dims for current curve
    // ie, we are finding weights of high-d points for a 2-d curve
    int pt_dim = mfa.domain.cols();             // dimensionality of input and control points (domain and range)
    vector<MatrixX<T>> NtQ(pt_dim);             // temp. matrices N^T x Q  for each dim of points (TODO: any way to avoid?)
    vector<MatrixX<T>> NtQN(pt_dim);            // matrices N^T x Q x N for each dim of points

    // allocate matrices of NtQ, NtQ2, NtQ2N, and NtQN
    for (auto k = 0; k < pt_dim; k++)
    {
        NtQ[k]   = MatrixX<T>::Zero(Nt.rows(),  Nt.cols());     // rectangular, n - 1 x m - 1
        NtQN[k]  = MatrixX<T>::Zero(N.rows(), N.rows());        // square, n -1 x n - 1
    }

    // compute NtQN and the columns of R for each dimension
    for (auto k = 0; k < pt_dim; k++)           // for all point dims
    {
        for (auto i = 0; i < Nt.cols(); i++)
        {
            T dom_pt_coord = mfa.domain(co + (i + 1) * mfa.ds[cur_dim], k); // current coord. of current input pt.
            NtQ[k].col(i)  = Nt.col(i) * dom_pt_coord;
        }
        NtQN[k]  = NtQ[k] * N;

        // compute R
        R.col(k) = NtQN[k] * weights.segment(1, N.cols());
    }

    // debug
    cerr << "R:\n" << R << endl;
}

// computes weighted right hand side vector of M&K'95 eq. 5 premultiplied by N^T, ie R = N^T * Q * N * w
// input is a new set of input points, not the default domain
//
// where:
// N is the matrix of basis function coefficients
// Q is a diagonal matrix of input points on the diagonals and 0 elsewhere (excluding first and last points)
// w is the weights vector, excluding first and last weights
// R is column vector of n - 1 elements, each element multiple coordinates of the input points
// ie, a matrix of n - 1 rows and control point dims columns
template <typename T>
void
mfa::
Encoder<T>::
RHS(
        int         cur_dim,             // current dimension
        MatrixX<T>& in_pts,              // input points (not the default domain stored in the mfa)
        MatrixX<T>& N,                   // matrix of basis function coefficients
        MatrixX<T>& Nt,                  // transpose of N
        MatrixX<T>& R,                   // (output) residual matrix allocated by caller
        VectorX<T>& weights,             // precomputed weights for n + 1 control points on this curve
        int         ko,                  // index of starting knot
        int         po,                  // index of starting parameter
        int         co,                  // index of starting domain pt in current curve
        int         cs)                  // stride of input pts in current curve
{
    // TODO: decide if I need to represent all dims or only 2 dims for current curve
    // ie, we are finding weights of high-d points for a 2-d curve
    int pt_dim = mfa.domain.cols();             // dimensionality of input and control points (domain and range)
    vector<MatrixX<T>> NtQ(pt_dim);             // temp. matrices N^T x Q  for each dim of points (TODO: any way to avoid?)
    vector<MatrixX<T>> NtQN(pt_dim);            // matrices N^T x Q x N for each dim of points

    // allocate matrices of NtQ, NtQ2, NtQ2N, and NtQN
    for (auto k = 0; k < pt_dim; k++)
    {
        NtQ[k]   = MatrixX<T>::Zero(Nt.rows(),  Nt.cols());     // rectangular, n - 1 x m - 1
        NtQN[k]  = MatrixX<T>::Zero(N.rows(), N.rows());        // square, n -1 x n - 1
    }

    // compute NtQN and the columns of R for each dimension
    for (auto k = 0; k < pt_dim; k++)           // for all point dims
    {
        for (auto i = 0; i < Nt.cols(); i++)
        {
            T dom_pt_coord = in_pts(co + (i + 1) * cs, k);   // current coord. of current input pt.
            NtQ[k].col(i)  = Nt.col(i) * dom_pt_coord;
        }
        NtQN[k]  = NtQ[k] * N;

        // compute R
        R.col(k) = NtQN[k] * weights.segment(1, N.cols());
    }
}

// Checks quantities needed for approximation
template <typename T>
void
mfa::
Encoder<T>::
Quants(
        VectorXi& n,                // (output) number of control point spans in each dim
        VectorXi& m)                // (output) number of input data point spans in each dim
{
    if (mfa.p.size() != mfa.ndom_pts.size())
    {
        fprintf(stderr, "Error: Encode() size of p must equal size of ndom_pts\n");
        exit(1);
    }
    for (size_t i = 0; i < mfa.p.size(); i++)
    {
        if (mfa.nctrl_pts(i) <= mfa.p(i))
        {
            fprintf(stderr, "Error: Encode() number of control points in dimension %ld"
                    "must be at least p + 1 for dimension %ld\n", i, i);
            exit(1);
        }
        if (mfa.nctrl_pts(i) > mfa.ndom_pts(i))
        {
            fprintf(stderr, "Warning: Encode() number of control points (%d) in dimension %ld "
                    "exceeds number of input data points (%d) in dimension %ld. "
                    "Technically, this is not an error, but it could be a sign something is wrong and "
                    "probably not desired if you want compression. You may not be able to get the "
                    "desired error limit and compression simultaneously. Try increasing error limit?\n",
                    mfa.nctrl_pts(i), i, mfa.ndom_pts(i), i);
//             exit(1);
        }
    }

    n.resize(mfa.p.size());
    m.resize(mfa.p.size());
    for (size_t i = 0; i < mfa.p.size(); i++)
    {
        n(i)        =  mfa.nctrl_pts(i) - 1;
        m(i)        =  mfa.ndom_pts(i)  - 1;
    }
}

// append points from P to temporary control points
// init first and last control points and copy rest from solution P
// TODO: any way to avoid this copy?
// last dimension gets copied to final control points
// previous dimensions get copied to alternating double buffers
template <typename T>
void
mfa::
Encoder<T>::
CopyCtrl(MatrixX<T>& P,          // solved points for current dimension and curve
         int         k,          // current dimension
         size_t      co,         // starting offset for reading domain points
         size_t      cs,         // stride for reading domain points
         size_t      to,         // starting offset for writing control points
         MatrixX<T>& temp_ctrl0, // first temporary control points buffer
         MatrixX<T>& temp_ctrl1) // second temporary control points buffer
{
    int ndims = mfa.ndom_pts.size();             // number of domain dimensions

    // if there is only one dim, copy straight to output
    if (ndims == 1)
    {
        mfa.ctrl_pts.row(to) = mfa.domain.row(co);
        for (int i = 1; i < P.rows() + 1; i++)
            mfa.ctrl_pts.row(to + i * cs) = P.row(i - 1);
        mfa.ctrl_pts.row(to + (P.rows() + 1) * cs) = mfa.domain.row(co + mfa.ndom_pts(k) - 1);
    }
    // first dim copied from domain to temp_ctrl0
    else if (k == 0)
    {
        temp_ctrl0.row(to) = mfa.domain.row(co);
        for (int i = 1; i < P.rows() + 1; i++)
            temp_ctrl0.row(to + i * cs) = P.row(i - 1);
        temp_ctrl0.row(to + (P.rows() + 1) * cs) = mfa.domain.row(co + mfa.ndom_pts(k) - 1);
    }
    // even numbered dims (but not the last one) copied from temp_ctrl1 to temp_ctrl0
    else if (k % 2 == 0 && k < ndims - 1)
    {
        temp_ctrl0.row(to) = temp_ctrl1.row(co);
        for (int i = 1; i < P.rows() + 1; i++)
            temp_ctrl0.row(to + i * cs) = P.row(i - 1);
        temp_ctrl0.row(to + (P.rows() + 1) * cs) = temp_ctrl1.row(co + (mfa.ndom_pts(k) - 1) * cs);
    }
    // odd numbered dims (but not the last one) copied from temp_ctrl0 to temp_ctrl1
    else if (k % 2 == 1 && k < ndims - 1)
    {
        temp_ctrl1.row(to) = temp_ctrl0.row(co);
        for (int i = 1; i < P.rows() + 1; i++)
            temp_ctrl1.row(to + i * cs) = P.row(i - 1);
        temp_ctrl1.row(to + (P.rows() + 1) * cs) = temp_ctrl0.row(co + (mfa.ndom_pts(k) - 1) * cs);
    }
    // final dim if even is copied from temp_ctrl1 to ctrl_pts
    else if (k == ndims - 1 && k % 2 == 0)
    {
        mfa.ctrl_pts.row(to) = temp_ctrl1.row(co);
        for (int i = 1; i < P.rows() + 1; i++)
            mfa.ctrl_pts.row(to + i * cs) = P.row(i - 1);
        mfa.ctrl_pts.row(to + (P.rows() + 1) * cs) = temp_ctrl1.row(co + (mfa.ndom_pts(k) - 1) * cs);
    }
    // final dim if odd is copied from temp_ctrl0 to ctrl_pts
    else if (k == ndims - 1 && k % 2 == 1)
    {
        mfa.ctrl_pts.row(to) = temp_ctrl0.row(co);
        for (int i = 1; i < P.rows() + 1; i++)
            mfa.ctrl_pts.row(to + i * cs) = P.row(i - 1);
        mfa.ctrl_pts.row(to + (P.rows() + 1) * cs) = temp_ctrl0.row(co + (mfa.ndom_pts(k) - 1) * cs);
    }
}

// append points from P to temporary control points
// init first and last control points and copy rest from solution P
// TODO: any way to avoid this copy?
// just simple copy to one temporary buffer, no alternating double buffers
// nor copy to final control points
template <typename T>
void
mfa::
Encoder<T>::
CopyCtrl(MatrixX<T>& P,          // solved points for current dimension and curve
         int         k,          // current dimension
         size_t      co,         // starting offset for reading domain points
         MatrixX<T>& temp_ctrl)  // temporary control points buffer
{
    // copy first point straight from domain
    temp_ctrl.row(0) = mfa.domain.row(co);

    // copy intermediate points
    // clamp all dimensions other than k to the same as the domain points
    // this eliminates any wiggles in other dimensions from the computation of P (typically ~10^-5)
    for (int i = 1; i < P.rows() + 1; i++)
    {
        for (auto j = 0; j < mfa.domain.cols(); j++)
        {
            if (j < mfa.p.size() && j != k)
                temp_ctrl(i, j) = mfa.domain(co, j);
            else
                temp_ctrl(i, j) = P(i - 1, j);
        }
    }

    // copy last point straight from domain
    temp_ctrl.row(P.rows() + 1) = mfa.domain.row(co + (mfa.ndom_pts(k) - 1) * mfa.ds[k]);
}

// solves for one curve of control points
template <typename T>
void
mfa::
Encoder<T>::
CtrlCurve(MatrixX<T>& N,          // basis functions for current dimension
          MatrixX<T>& Nt,         // transpose of N
          MatrixX<T>& NtN,        // Nt * N
          MatrixX<T>& NtNi,       // inverse of NtN
          MatrixX<T>& R,          // residual matrix for current dimension and curve
          MatrixX<T>& P,          // solved points for current dimension and curve
          size_t      k,          // current dimension
          size_t      co,         // starting ofst for reading domain pts
          size_t      cs,         // stride for reading domain points
          size_t      to,         // starting ofst for writing control pts
          MatrixX<T>& temp_ctrl0, // first temporary control points buffer
          MatrixX<T>& temp_ctrl1) // second temporary control points buffer
{
    int n = N.cols() + 1;               // number of control point spans
    int m = N.rows() + 1;               // number of input point spans

    // solve for weights
    // TODO: avoid copying into Q by passing temp_ctrl0, temp_ctrl1, co, cs to Weights()
    // TODO: check that this is right, using co and cs for copying control points and domain points
    MatrixX<T> Q;
    Q.resize(mfa.ndom_pts(k), mfa.domain.cols());
    if (k == 0)
    {
        for (auto i = 0; i < mfa.ndom_pts(k); i++)
            Q.row(i) = mfa.domain.row(co + i * cs);
    }
    else if (k % 2)
    {
        for (auto i = 0; i < mfa.ndom_pts(k); i++)
            Q.row(i) = temp_ctrl0.row(co + i * cs);
    }
    else
    {
        for (auto i = 0; i < mfa.ndom_pts(k); i++)
            Q.row(i) = temp_ctrl1.row(co + i * cs);
    }

    // debug
//     cerr << "k=" << k << " Q:\n" << Q << endl;

    // solve for weights
    VectorX<T> weights = VectorX<T>::Ones(n + 1);    // weights for control points on this curve
    Weights(Q, N, Nt, NtN, NtNi, weights);

    // compute R
    // first dimension reads from domain
    // subsequent dims alternate reading temp_ctrl0 and temp_ctrl1
    // even dim reads temp_ctrl1, odd dim reads temp_ctrl0; opposite of writing order
    // because what was written in the previous dimension is read in the current one

#if 1       // RHS from P&T

    if (k == 0)
        RHS(k, N, R, weights, mfa.ko[k], mfa.po[k], co);                 // input points = default domain
    else if (k % 2)
        RHS(k, temp_ctrl0, N, R, weights, mfa.ko[k], mfa.po[k], co, cs); // input points = temp_ctrl0
    else
        RHS(k, temp_ctrl1, N, R, weights, mfa.ko[k], mfa.po[k], co, cs); // input points = temp_ctrl1

#else       // RHS from M&K

    if (k == 0)
        RHS(k, N, Nt, R, weights, mfa.ko[k], mfa.po[k], co);                 // input points = default domain
    else if (k % 2)
        RHS(k, temp_ctrl0, N, Nt,  R, weights, mfa.ko[k], mfa.po[k], co, cs); // input points = temp_ctrl0
    else
        RHS(k, temp_ctrl1, N, Nt,  R, weights, mfa.ko[k], mfa.po[k], co, cs); // input points = temp_ctrl1

#endif

//     // compute rational denominators for input points
//     // rational denominator requires all n + 1 basis functions and weights, not skipping first and last
//     VectorX<T> denom(m - 1);            // rational denomoninator for param of each input point excl. ends
//     MatrixX<T> Nk(1, n + 1);            // basis function coeffs for one parameter value (all coeffs, not skipping first and last)
//     for (int j = 1; j < m; j++)
//     {
//         int span = mfa.FindSpan(k, mfa.params(mfa.po[k] + j), mfa.ko[k]) - mfa.ko[k];   // relative to ko
//         assert(span <= n);                   // sanity
//         Nk = MatrixX<T>::Zero(1, n + 1);     // basis coefficients for Rk[j]
//         mfa.BasisFuns(k, mfa.params(mfa.po[k] + j), span, Nk, 0, n, 0);
//         denom(j - 1) = (Nk.row(0).cwiseProduct(weights.transpose())).sum();
//     }
// 
//     //         cerr << "denom:\n" << denom << endl;
// 
//     // "rationalize" N and Nt
//     // ie, convert their basis function coefficients to rational ones with weights
//     MatrixX<T> N_rat = N;                       // need a copy because N, Nt will be reused for other curves
//     MatrixX<T> Nt_rat = Nt;
//     MatrixX<T> NtN_rat = NtN;
//     for (auto i = 0; i < N.cols(); i++)
//         N_rat.col(i) *= weights(i + 1);
//     for (auto j = 0; j < N.rows(); j++)
//         N_rat.row(j) /= denom(j);
//     for (auto i = 0; i < Nt.rows(); i++)
//         Nt_rat.row(i) *= weights(i + 1);
//     for (auto j = 0; j < Nt.cols(); j++)
//         Nt_rat.col(j) /= denom(j);
//     NtN_rat = Nt_rat * N_rat;

    // debug
    //         cerr << "k " << k << " NtN:\n" << NtN << endl;
    //         cerr << " NtN_rat:\n" << NtN_rat << endl;

    // rationalize NtN, ie, weight the basis function coefficients
    MatrixX<T> NtN_rat = NtN;
    mfa.Rationalize(k, weights, N, Nt, NtN_rat);

    // solve for P
    P = NtN_rat.ldlt().solve(R);

    // append points from P to control points
    // TODO: any way to avoid this?
    CopyCtrl(P, k, co, cs, to, temp_ctrl0, temp_ctrl1);
}

// returns number of points in a curve that have error greater than err_limit
template <typename T>
int
mfa::
Encoder<T>::
ErrorCurve(
        size_t         k,                         // current dimension
        size_t         co,                        // starting ofst for reading domain pts
        MatrixX<T>&    ctrl_pts,                  // control points
        VectorX<T>&    weights,                   // weights associated with control points
        T          err_limit)                 // max allowable error
{
    mfa::Decoder<T> decoder(mfa);
    VectorX<T> cpt(mfa.domain.cols());            // decoded curve point
    int nerr = 0;                               // number of points with error greater than err_limit
    int span = mfa.p[k];                        // current knot span of the domain point being checked

    for (auto i = 0; i < mfa.ndom_pts[k]; i++)      // all domain points in the curve
    {
        while (mfa.knots(mfa.ko[k] + span + 1) < 1.0 && mfa.knots(mfa.ko[k] + span + 1) <= mfa.params(mfa.po[k] + i))
            span++;
        // debug
//         fprintf(stderr, "param=%.3f span=[%.3f %.3f]\n", mfa.params(po[k] + i), knots(mfa.ko[k] + span), knots(mfa.ko[k] + span + 1));

        decoder.CurvePt(k, mfa.params(mfa.po[k] + i), ctrl_pts, weights, cpt, mfa.ko[k]);
        T err = fabs(mfa.NormalDistance(cpt, co + i * mfa.ds[k])) / mfa.range_extent;       // normalized by data range
//         T err = fabs(mfa.CurveDistance(k, cpt, co + i * mfa.ds[k])) / mfa.dom_range;     // normalized by data range
        if (err > err_limit)
        {
            nerr++;

            // debug
//             VectorX<T> dpt = mfa.domain.row(co + i * mfa.ds[k]);
//             cerr << "\ndomain point:\n" << dpt << endl;
//             cerr << "approx point:\n" << cpt << endl;
//             fprintf(stderr, "k=%ld i=%d co=%ld err=%.3e\n\n", k, i, co, err);
        }
    }

    return nerr;
}

// computes new knots to be inserted into a curve
// for each current knot span where the error is greater than the limit, finds the domain point
// where the error is greatest and adds the knot at that parameter value
//
// this version takes a set of control points as input instead of mfa.ctrl_pts
template <typename T>
void
mfa::
Encoder<T>::
ErrorCurve(
        size_t           k,                       // current dimension
        size_t           co,                      // starting ofst for reading domain pts
        MatrixX<T>&      ctrl_pts,                // control points
        VectorX<T>&      weights,                 // weights associated with control points
        VectorXi&        nnew_knots,              // number of new knots
        vector<T>&       new_knots,               // new knots
        T                err_limit)               // max allowable error
{
    mfa::Decoder<T> decoder(mfa);
    VectorX<T> cpt(mfa.domain.cols());            // decoded curve point
    int span      = mfa.p[k];                    // current knot span of the domain point being checked
    int old_span  = -1;                          // span of previous domain point
    T max_err = 0;                          // max error seen so far in the same span
    size_t max_err_pt;                          // index of domain point in same span with max error
    bool new_split = false;                     // a new split was found in the current span

    for (auto i = 0; i < mfa.ndom_pts[k]; i++)      // all domain points in the curve
    {
        while (mfa.knots(mfa.ko[k] + span + 1) < 1.0 && mfa.knots(mfa.ko[k] + span + 1) <= mfa.params(mfa.po[k] + i))
            span++;

        if (span != old_span)
            max_err = 0;

        // record max of previous span if span changed and previous span had a new split
        if (span != old_span && new_split)
        {
            nnew_knots(k)++;
            new_knots.push_back(mfa.params(mfa.po[k] + max_err_pt));
            new_split = false;
        }

        decoder.CurvePt(k, mfa.params(mfa.po[k] + i), ctrl_pts, weights, cpt, mfa.ko[k]);

        T err = fabs(mfa.NormalDistance(cpt, co + i * mfa.ds[k])) / mfa.range_extent;     // normalized by data range

        if (err > err_limit && err > max_err)  // potential new knot
        {
            // ensure there would be a domain point in both halves of the span if it were split
            bool split_left = false;
            for (auto j = i; mfa.params(mfa.po[k] + j) >= mfa.knots(mfa.ko[k] + span); j--)
                if (mfa.params(mfa.po[k] + j) < mfa.params(mfa.po[k] + i))
                {
                    split_left = true;
                    break;
                }
            bool split_right = false;
            for (auto j = i; mfa.params(mfa.po[k] + j) < mfa.knots(mfa.ko[k] + span + 1); j++)
                if (mfa.params(mfa.po[k] + j) >= mfa.params(mfa.po[k] + i))
                {
                    split_right = true;
                    break;
                }
            // record the potential split point
            if (split_left && split_right && err > max_err)
            {
                max_err = err;
                max_err_pt = i;
                new_split = true;
            }
        }                                                           // potential new knot

        if (span != old_span)
            old_span = span;
    }

    // record max of last span
    if (new_split)
    {
        nnew_knots(k)++;
        new_knots.push_back(mfa.params(mfa.po[k] + max_err_pt));
    }
}

// returns number of points in a curve that have error greater than err_limit
// fills err_spans with the span indices of spans that have at least one point with such error
//  and that have at least one inut point in each half of the span (assuming eventually
//  the span would be split in half with a knot added in the middle, and an input point would
//  need to be in each span after splitting)
//
// this version takes a set instead of a vector for error_spans so that the same span can be
// added iteratively multiple times without creating duplicates
//
// this version takes a set of control points as input instead of mfa.ctrl_pts
template <typename T>
int
mfa::
Encoder<T>::
ErrorCurve(
        size_t         k,                         // current dimension
        size_t         co,                        // starting ofst for reading domain pts
        MatrixX<T>&    ctrl_pts,                  // control points
        VectorX<T>&    weights,                   // weights associated with control points
        set<int>&      err_spans,                 // spans with error greater than err_limit
        T              err_limit)                 // max allowable error
{
    mfa::Decoder<T> decoder(mfa);
    VectorX<T> cpt(mfa.domain.cols());            // decoded curve point
    int nerr = 0;                               // number of points with error greater than err_limit
    int span = mfa.p[k];                        // current knot span of the domain point being checked

    for (auto i = 0; i < mfa.ndom_pts[k]; i++)      // all domain points in the curve
    {
        while (mfa.knots(mfa.ko[k] + span + 1) < 1.0 && mfa.knots(mfa.ko[k] + span + 1) <= mfa.params(mfa.po[k] + i))
            span++;

        decoder.CurvePt(k, mfa.params(mfa.po[k] + i), ctrl_pts, weights, cpt, mfa.ko[k]);

        T err = fabs(mfa.NormalDistance(cpt, co + i * mfa.ds[k])) / mfa.range_extent;       // normalized by data range
//         T err = fabs(mfa.CurveDistance(k, cpt, co + i * mfa.ds[k])) / mfa.dom_range;     // normalized by data range

        if (err > err_limit)
        {
            // don't duplicate spans
            set<int>::iterator it = err_spans.find(span);
            if (!err_spans.size() || it == err_spans.end())
            {
                // ensure there would be a domain point in both halves of the span if it were split
                bool split_left = false;
                for (auto j = i; mfa.params(mfa.po[k] + j) >= mfa.knots(mfa.ko[k] + span); j--)
                    if (mfa.params(mfa.po[k] + j) < (mfa.knots(mfa.ko[k] + span) + mfa.knots(mfa.ko[k] + span + 1)) / 2.0)
                    {
                        // debug
//                         fprintf(stderr, "split_left: param=%.3f span[%d]=[%.3f, %.3f)\n",
//                                 mfa.params(po[k] + j), span, knots(mfa.ko[k] + span), knots(mfa.ko[k] + span + 1));
                        split_left = true;
                        break;
                    }
                bool split_right = false;
                for (auto j = i; mfa.params(mfa.po[k] + j) < mfa.knots(mfa.ko[k] + span + 1); j++)
                    if (mfa.params(mfa.po[k] + j) >= (mfa.knots(mfa.ko[k] + span) + mfa.knots(mfa.ko[k] + span + 1)) / 2.0)
                    {
//                         fprintf(stderr, "split_right: param=%.3f span[%d]=[%.3f, %.3f)\n",
//                                 mfa.params(mfa.po[k] + j), span, knots(mfa.ko[k] + span), knots(mfa.ko[k] + span + 1));
                        split_right = true;
                        break;
                    }
                // mark the span and count the point if the span can (later) be split
                if (split_left && split_right)
                    err_spans.insert(it, span);
            }
            // count the point in the total even if the span is not marked for splitting
            // total used to find worst curve, defined as the curve with the most domain points in
            // error (including multiple domain points per span and points in spans that can't be
            // split further)
            nerr++;

            // debug
//             VectorX<T> dpt = mfa.domain.row(co + i * mfa.ds[k]);
//             cerr << "\ndomain point:\n" << dpt << endl;
//             cerr << "approx point:\n" << cpt << endl;
//             fprintf(stderr, "k=%ld i=%d co=%ld err=%.3e\n\n", k, i, co, err);
        }
    }

    return nerr;
}

// returns number of points in a curve that have error greater than err_limit
// fills err_spans with the span indices of spans that have at least one point with such error
//  and that have at least one inut point in each half of the span (assuming eventually
//  the span would be split in half with a knot added in the middle, and an input point would
//  need to be in each span after splitting)
//
// this version takes a set instead of a vector for error_spans so that the same span can be
// added iteratively multiple times without creating duplicates
//
// this version uses mfa.ctrl_pts for control points
template <typename T>
int
mfa::
Encoder<T>::
ErrorCurve(
        size_t       k,                         // current dimension
        size_t       co,                        // starting ofst for reading domain pts
        size_t       to,                        // starting ofst for reading control pts
        set<int>&    err_spans,                 // spans with error greater than err_limit
        T            err_limit)                 // max allowable error
{
    mfa::Decoder<T> decoder(mfa);
    VectorX<T> cpt(mfa.domain.cols());            // decoded curve point
    int nerr = 0;                               // number of points with error greater than err_limit
    int span = mfa.p[k];                        // current knot span of the domain point being checked

    for (auto i = 0; i < mfa.ndom_pts[k]; i++)      // all domain points in the curve
    {
        while (mfa.knots(mfa.ko[k] + span + 1) < 1.0 && mfa.knots(mfa.ko[k] + span + 1) <= mfa.params(mfa.po[k] + i))
            span++;

        decoder.CurvePt(k, mfa.params(mfa.po[k] + i), to, cpt);

        T err = fabs(mfa.NormalDistance(cpt, co + i * mfa.ds[k])) / mfa.range_extent;     // normalized by data range

        if (err > err_limit)
        {
            // don't duplicate spans
            set<int>::iterator it = err_spans.find(span);
            if (!err_spans.size() || it == err_spans.end())
            {
                // ensure there would be a domain point in both halves of the span if it were split
                bool split_left = false;
                for (auto j = i; mfa.params(mfa.po[k] + j) >= mfa.knots(mfa.ko[k] + span); j--)
                    if (mfa.params(mfa.po[k] + j) < (mfa.knots(mfa.ko[k] + span) + mfa.knots(mfa.ko[k] + span + 1)) / 2.0)
                    {
                        split_left = true;
                        break;
                    }
                bool split_right = false;
                for (auto j = i; mfa.params(mfa.po[k] + j) < mfa.knots(mfa.ko[k] + span + 1); j++)
                    if (mfa.params(mfa.po[k] + j) >= (mfa.knots(mfa.ko[k] + span) + mfa.knots(mfa.ko[k] + span + 1)) / 2.0)
                    {
                        split_right = true;
                        break;
                    }
                // mark the span and count the point if the span can (later) be split
                if (split_left && split_right)
                    err_spans.insert(it, span);
            }
            // count the point in the total even if the span is not marked for splitting
            nerr++;
        }
    }

    return nerr;
}

// error of points decoded from a curve aligned with a curve of control points
//
// returns number of points in a curve that have error greater than err_limit
// fills err_spans with the span indices of spans that have at least one point with such error
//  and that have at least one inut point in each half of the span (assuming eventually
//  the span would be split in half with a knot added in the middle, and an input point would
//  need to be in each span after splitting)
//
// this version takes a set instead of a vector for error_spans so that the same span can be
// added iteratively multiple times without creating duplicates
//
// this version uses mfa.ctrl_pts for control points
template <typename T>
int
mfa::
Encoder<T>::
ErrorCtrlCurve(
        size_t       k,                         // current dimension
        size_t       to,                        // starting ofst for reading control pts
        set<int>&    err_spans,                 // spans with error greater than err_limit
        T            err_limit)                 // max allowable error
{
    mfa::Decoder<T> decoder(mfa);
    VectorX<T> cpt(mfa.domain.cols());            // decoded curve point
    int nerr = 0;                               // number of points with error greater than err_limit
    int span = mfa.p[k];                        // current knot span of the domain point being checked

    // compute parameter value of start of control curve
    vector<T> param(mfa.p.size());
    for (auto k = 0; k < mfa.p.size(); k++)
        // TODO: decide whether to use InterpolateParams() or write Param() function for one target point
        // (InterpolateParams has assumption of increasing domain coordinates that needs to be removed)
        param[k] = mfa.InterpolateParams(k, mfa.po[k], mfa.ds[k], mfa.ctrl_pts(to, k));

    // debug
//     fprintf(stderr, "param = [ ");
//     for (auto k = 0; k < mfa.p.size(); k++)
//         fprintf(stderr, "%.3f ", param[k]);
//     fprintf(stderr, "]\n");

    // compute ijk index of closest params to parameter value at start of control curve
    VectorXi ijk(mfa.p.size());
    for (auto k = 0; k < mfa.p.size(); k++)
    {
        size_t j;
        for (j = 0; j < mfa.ndom_pts[k] && mfa.params[mfa.po[k] + j] < param[k]; j++)
            ;
        ijk(k) = (j > 0 ? j - 1 : j);
    }

    size_t co;              // starting offset of domain points for curve closest to control point curve
    mfa.ijk2idx(ijk, co);

    // debug
//     VectorX<T> ctpt = mfa.ctrl_pts.row(to);
//     VectorX<T> dopt = mfa.domain.row(co);
//     cerr << "start ctrl pt:\n" << ctpt << endl;
//     cerr << "start inpt pt:\n" << dopt << "\n" << endl;

    for (auto i = 0; i < mfa.ndom_pts[k]; i++)      // all domain points in the curve
    {
        while (mfa.knots(mfa.ko[k] + span + 1) < 1.0 && mfa.knots(mfa.ko[k] + span + 1) <= mfa.params(mfa.po[k] + i))
            span++;

        decoder.CurvePt(k, mfa.params(mfa.po[k] + i), to, cpt);

        // adjust input point index so that input point cell contains cpt
        size_t j = i;
        while (j < mfa.ndom_pts[k] - 1 && mfa.domain(co + j * mfa.ds[k], k) < cpt(k))
            j++;
        while (j > 0 &&                   mfa.domain(co + j * mfa.ds[k], k) > cpt(k))
            j--;

        // debug: check that the input cell contains cpt
        bool error = false;
        if (j == 0 && mfa.domain(co + (j + 1) * mfa.ds[k], k) < cpt(k))
            error = true;
        else if (j == mfa.ndom_pts[k] - 2 && mfa.domain(co + j * mfa.ds[k], k) > cpt(k))
            error = true;
        else if (j > 0 && j < mfa.ndom_pts[k] - 2 &&
                (mfa.domain(co + j * mfa.ds[k], k) > cpt(k) || mfa.domain(co + (j + 1) * mfa.ds[k], k) < cpt(k)))
            error = true;
        if (error)
            fprintf(stderr, "Error: j=%ld %f is not contained in [%f, %f]\n", j,
                    cpt(k), mfa.domain(co + j * mfa.ds[k], k), mfa.domain(co + (j + 1) * mfa.ds[k], k));

        // compute error
//         T err = fabs(mfa.NormalDistance(cpt, co + j * mfa.ds[k])) / mfa.range_extent;     // normalized by data range
        T err = fabs(mfa.CurveDistance(k, cpt, co + j * mfa.ds[k])) / mfa.range_extent;      // normalized by data range

        if (err > err_limit)
        {
            // debug
//             VectorX<T> dopt = mfa.domain.row(co + j * mfa.ds[k]);
//             cerr << "decoded pt:\n" << cpt  << endl;
//             cerr << "input pt:\n"   << dopt << endl;
//             fprintf(stderr, "err = %.3e\n\n", err);

            // don't duplicate spans
            set<int>::iterator it = err_spans.find(span);
            if (!err_spans.size() || it == err_spans.end())
            {
                // ensure there would be a domain point in both halves of the span if it were split
                bool split_left = false;
                for (auto j = i; mfa.params(mfa.po[k] + j) >= mfa.knots(mfa.ko[k] + span); j--)
                    if (mfa.params(mfa.po[k] + j) < (mfa.knots(mfa.ko[k] + span) + mfa.knots(mfa.ko[k] + span + 1)) / 2.0)
                    {
                        split_left = true;
                        break;
                    }
                bool split_right = false;
                for (auto j = i; mfa.params(mfa.po[k] + j) < mfa.knots(mfa.ko[k] + span + 1); j++)
                    if (mfa.params(mfa.po[k] + j) >= (mfa.knots(mfa.ko[k] + span) + mfa.knots(mfa.ko[k] + span + 1)) / 2.0)
                    {
                        split_right = true;
                        break;
                    }
                // mark the span and count the point if the span can (later) be split
                if (split_left && split_right)
                    err_spans.insert(it, span);
            }
            // count the point in the total even if the span is not marked for splitting
            nerr++;
        }
    }

    return nerr;
}

#include    "encode_templates.cpp"
