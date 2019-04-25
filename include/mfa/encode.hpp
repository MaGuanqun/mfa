//--------------------------------------------------------------
// encoder object
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#ifndef _ENCODE_HPP
#define _ENCODE_HPP

#include    <mfa/data_model.hpp>
#include    <mfa/mfa.hpp>
#include    <mfa/decode.hpp>
#include    <mfa/new_knots.hpp>

#include    <Eigen/Dense>
#include    <vector>
#include    <set>
#include    <iostream>

#ifndef      MFA_NO_WEIGHTS

#include    "coin/ClpSimplex.hpp"
#include    "coin/ClpInterior.hpp"

#endif

typedef Eigen::MatrixXf MatrixXf;
typedef Eigen::MatrixXi MatrixXi;
typedef Eigen::VectorXf VectorXf;
typedef Eigen::VectorXi VectorXi;

using namespace std;

template <typename T>                       // float or double
class NewKnots;


namespace mfa
{
    template <typename T>                   // float or double
    class Encoder
    {
    public:

        Encoder(
                MatrixX<T>&     domain_,            // input points
                MFA_Data<T>&    mfa_,               // MFA data model
                int             verbose_)           // output level
            : mfa(mfa_),
            verbose(verbose_),
            domain(domain_),
            max_num_curves(1.0e4)                           // max num. curves to check in one dimension of curve version
        {}

        ~Encoder() {}

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
        void Encode(bool weighted = true)                   // solve for and use weights
        {
            // check and assign main quantities
            VectorXi n;                                     // number of control point spans in each domain dim
            VectorXi m;                                     // number of input data point spans in each domain dim
            int      ndims      = mfa.ndom_pts.size();      // number of domain dimensions
            size_t   cs         = 1;                        // stride for input points in curve in cur. dim
            int      pt_dim     = mfa.tmesh.tensor_prods[0].ctrl_pts.cols();    // control point dimensonality

            Quants(n, m);

            // control points
            // TODO: hard-coded for single tensor
            mfa.tmesh.tensor_prods[0].ctrl_pts.resize(mfa.tmesh.tensor_prods[0].nctrl_pts.prod(), mfa.max_dim - mfa.min_dim + 1);

            // 2 buffers of temporary control points
            // double buffer needed to write output curves of current dim without changing its input pts
            // temporary control points need to begin with size as many as the input domain points
            // except for the first dimension, which can be the correct number of control points
            // because the input domain points are converted to control points one dimension at a time
            // TODO: need to find a more space-efficient way
            size_t tot_ntemp_ctrl = 1;
            for (size_t k = 0; k < ndims; k++)
            {
                // TODO: hard-coded for single tensor
                tot_ntemp_ctrl *= (k == 0 ? mfa.tmesh.tensor_prods[0].nctrl_pts(k) : mfa.ndom_pts(k));
            }
            MatrixX<T> temp_ctrl0 = MatrixX<T>::Zero(tot_ntemp_ctrl, pt_dim);
            MatrixX<T> temp_ctrl1 = MatrixX<T>::Zero(tot_ntemp_ctrl, pt_dim);

            VectorXi ntemp_ctrl = mfa.ndom_pts;     // current num of temp control pts in each dim

            for (size_t k = 0; k < ndims; k++)      // for all domain dimensions
            {
                // number of curves in this dimension
                size_t ncurves;
                ncurves = 1;
                for (int i = 0; i < ndims; i++)
                {
                    if (i < k)
                    {
                        // TODO: hard-coded for single tensor
                        ncurves *= mfa.tmesh.tensor_prods[0].nctrl_pts(i);
                    }
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
                        // TODO: hard-coded for single tensor
                        to[j] = too + cs * mfa.tmesh.tensor_prods[0].nctrl_pts(k);
                        too   = to[j];
                    }
                }

                // TODO:
                // Investigate whether in later dimensions, when input data points are replaced by
                // control points, need new knots and params computed.
                // In the next dimension, the coordinates of the dimension didn't change,
                // but the chord length did, as control points moved away from the data points in
                // the prior dim. Need to determine how much difference it would make to recompute
                // params and knots for the new input points
                // (moot when using domain decomposition)

                // N is a matrix of (m + 1) x (n + 1) scalars that are the basis function coefficients
                //  _                          _
                // |  N_0(u[0])   ... N_n(u[0]) |
                // |     ...      ...      ...  |
                // |  N_0(u[m])   ... N_n(u[m]) |
                //  -                          -
                // TODO: N is going to be very sparse when it is large: switch to sparse representation
                // N has semibandwidth < p  nonzero entries across diagonal
                MatrixX<T> N = MatrixX<T>::Zero(m(k) + 1, n(k) + 1); // coefficients matrix

                for (int i = 0; i < N.rows(); i++)              // the rows of N
                {
                    // TODO: hard-coded for single tensor
                    int span = mfa.FindSpan(k, mfa.params(mfa.po[k] + i), mfa.tmesh.tensor_prods[0]);
                    mfa.BasisFuns(k, mfa.params(mfa.po[k] + i), span, N, i);
                }

                // TODO: NtN is going to be very sparse when it is large: switch to sparse representation
                // NtN has semibandwidth < p + 1 nonzero entries across diagonal
                MatrixX<T> NtN  = N.transpose() * N;

                // debug
//                 cerr << "N:\n" << N << endl;
//                 cerr << "NtN:\n" << NtN << endl;

#ifndef MFA_NO_TBB                                  // TBB version

                parallel_for (size_t(0), ncurves, [&] (size_t j)      // for all the curves in this dimension
                        {
                        // debug
                        // fprintf(stderr, "j=%ld curve\n", j);

                        // R is the right hand side needed for solving NtN * P = R
                        MatrixX<T> R(N.cols(), pt_dim);

                        // P are the unknown control points and the solution to NtN * P = R
                        // NtN is positive definite -> do not need pivoting
                        // TODO: use a common representation for P and ctrl_pts to avoid copying
                        MatrixX<T> P(N.cols(), pt_dim);

                        // compute the one curve of control points
                        // TODO: hard-coded for first tensor product
                        CtrlCurve(N, NtN, R, P, k, co[j], cs, to[j], temp_ctrl0, temp_ctrl1, -1, mfa.tmesh.tensor_prods[0], weighted);
                        });                                                  // curves in this dimension

#else                                       // serial vesion

                // R is the right hand side needed for solving NtN * P = R
                MatrixX<T> R(N.cols(), pt_dim);

                // P are the unknown control points and the solution to NtN * P = R
                // NtN is positive definite -> do not need pivoting
                // TODO: use a common representation for P and ctrl_pts to avoid copying
                MatrixX<T> P(N.cols(), pt_dim);

                // encode curves in this dimension
                for (size_t j = 0; j < ncurves; j++)
                {
                    // print progress
                    if (verbose)
                    {
                        if (j > 0 && j > 100 && j % (ncurves / 100) == 0)
                            fprintf(stderr, "\r dimension %ld: %.0f %% encoded (%ld out of %ld curves)",
                                    k, (T)j / (T)ncurves * 100, j, ncurves);
                    }

                    // compute the one curve of control points
                    // TODO: hard-coded for first tensor product
                    CtrlCurve(N, NtN, R, P, k, co[j], cs, to[j], temp_ctrl0, temp_ctrl1, j, mfa.tmesh.tensor_prods[0], weighted);
                }

#endif
                // adjust offsets and strides for next dimension
                // TODO: hard-coded for single tensor
                ntemp_ctrl(k) = mfa.tmesh.tensor_prods[0].nctrl_pts(k);
                cs *= ntemp_ctrl(k);

                NtN.resize(0, 0);                           // free NtN

                // print progress
                if (verbose)
                    fprintf(stderr, "\ndimension %ld of %d encoded\n", k + 1, ndims);

            }                                                      // domain dimensions
        }

        // adaptive encoding
        void AdaptiveEncode(
                T           err_limit,              // maximum allowable normalized error
                bool        weighted,               // solve for and use weights
                VectorX<T>& extents,                // extents in each dimension, for normalizing error (size 0 means do not normalize)
                int         max_rounds = 0)         // optional maximum number of rounds
        {
            VectorXi  nnew_knots = VectorXi::Zero(mfa.p.size());        // number of new knots in each dim
            vector<T> new_knots;                                        // new knots (1st dim changes fastest)

            // TODO: use weights for knot insertion
            // for now, weights are only used for final full encode

            // loop until no change in knots
            for (int iter = 0; ; iter++)
            {
                if (max_rounds > 0 && iter >= max_rounds)               // optional cap on number of rounds
                    break;

                if (verbose)
                    fprintf(stderr, "Iteration %d...\n", iter);

#ifdef HIGH_D           // high-d w/ splitting spans in the middle

                bool done = NewKnots_full(nnew_knots, new_knots, err_limit, extents, iter);

#else               // low-d w/ splitting spans in the middle

                bool done = NewKnots_curve(nnew_knots, new_knots, err_limit, extents, iter);

#endif

                // no new knots to be added
                if (done)
                {
                    if (verbose)
                        fprintf(stderr, "\nKnot insertion done after %d iterations; no new knots added.\n\n", iter + 1);
                    break;
                }

                // check if the new knots would make the number of control points >= number of input points in any dim
                done = false;
                for (auto k = 0; k < mfa.p.size(); k++)
                {
                    // TODO: hard-coded for single tensor
                    if (mfa.ndom_pts(k) <= mfa.tmesh.tensor_prods[0].nctrl_pts(k) + nnew_knots(k))
                    {
                        done = true;
                        break;
                    }
                }
                if (done)
                {
                    if (verbose)
                        fprintf(stderr, "\nKnot insertion done after %d iterations; control points would outnumber input points.\n", iter + 1);
                    break;
                }
            }

            // final full encoding needed after last knot insertion above
            if (verbose)
                fprintf(stderr, "Encoding in full %ldD\n", mfa.p.size());
            Encode(weighted);
        }

    private:

#ifndef      MFA_NO_WEIGHTS

        bool Weights(
                int         k,              // current dimension
                MatrixX<T>& Q,              // input points
                MatrixX<T>& N,              // basis functions
                MatrixX<T>& NtN,            // N^T * N
                int         curve_id,       // debugging
                VectorX<T>& weights)        // output weights
        {
            bool success;

            // Nt, NtNi
            // TODO: offer option of saving time or space by comuting Nt and NtN each time it is needed?
            MatrixX<T> Nt   = N.transpose();
            MatrixX<T> NtNi = NtN.partialPivLu().inverse();

            int pt_dim = mfa.max_dim - mfa.min_dim + 1;                     // dimensionality of input and control points (domain and range)

            MatrixX<T> NtQ2  = MatrixX<T>::Zero(Nt.rows(),   Nt.cols());    // N^T * Q^2
            MatrixX<T> NtQ   = MatrixX<T>::Zero(Nt.rows(),   Nt.cols());    // N^T * Q
            MatrixX<T> NtQ2N = MatrixX<T>::Zero(NtN.rows(),  NtN.cols());   // N^T * Q^2 * N
            MatrixX<T> NtQN  = MatrixX<T>::Zero(NtN.rows(),  NtN.cols());   // N^T * Q   * N

            // temporary matrices NtQ and NtQ2
            for (auto i = 0; i < Nt.cols(); i++)
            {
                NtQ.col(i)  = Nt.col(i) * Q(i, pt_dim - 1);
                NtQ2.col(i) = Nt.col(i) * Q(i, pt_dim - 1) * Q(i, pt_dim - 1);
            }

            // final matrices NtQN and NtQ2N
            NtQN  = NtQ  * N;
            NtQ2N = NtQ2 * N;

            // compute the matrix M according to eq.3 and eq. 4 of M&K95
            MatrixX<T> M = NtQ2N - NtQN * NtNi * NtQN;

            // debug: output the matrix M
            //     cerr << M << endl;
            //     Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
            //     ofstream M_out;
            //     M_out.open("M.txt");
            //     M_out << M.format(OctaveFmt) << endl;
            //     M_out.close();

            // compute the eigenvalues and eigenvectors of M (eq. 9 of M&K95)
            // NB: M is self-adjoint (ie, real and symmetric)
            // Eigen's SelfAdjointEigenSolver is faster than its normal EigenSolver in this case
            Eigen::SelfAdjointEigenSolver<MatrixX<T>> eigensolver(M);
            if (eigensolver.info() != Eigen::Success)
            {
                fprintf(stderr, "Error: Encoder::Weights(), computing eigenvalues of M failed, perhaps M is not self-adjoint?\n");
                return false;
            }

            const MatrixX<T>& EV    = eigensolver.eigenvectors();          // typing shortcut
            const VectorX<T>& evals = eigensolver.eigenvalues();           // typing shortcut

            // eigenvalues should be positive and distinct
            for (auto i = 0; i < evals.size() - 1; i++)
                if (evals(i) == 0.0 || evals(i) == evals(i + 1))
                {
                    fprintf(stderr, "Warning: Weights(): eigenvalues should be positive and distinct.\n");
                    fprintf(stderr, "Aborting weights calculation\n");
                    return false;
                }

            // if smallest eigenvector is all positive or all negative, those are the weights
            if ( (EV.col(0).array() > 0.0).all() )
            {
                weights = EV.col(0);
                weights *= (1.0 / weights.maxCoeff());  // scale to max weight = 1
                success = true;
            }
            else if ( (EV.col(0).array() < 0.0).all() )
            {
                weights = -EV.col(0);
                weights *= (1.0 / weights.maxCoeff());  // scale to max weight = 1
                success = true;
            }

            // if smallest eigenvector is mixed sign, then expand eigen space
            else
            {
                //         fprintf(stderr, "\nExpanding eigenspace using linear solver\n");
                success      = false;
                T min_weight = 1.0;
                T max_weight = 1.0e4;
                // minimum eigenvector element, if less, clamp to 0.0 (between 1e-6 and 1e-12 seems to help)
                T min_ev_val = 1.0e-12;

                auto nweights = weights.size();

                // column indices
                vector<int> col_idx(2 * nweights);
                for (auto j = 0; j < 2 * nweights; j++)
                    col_idx[j] = j;

                ClpSimplex model;                           // simplex method, ~10X faster than interior point
                //         ClpInterior model;                          // interior point method
                model.resize(2 * nweights, 0);              // all the rows and no columns (columns are added incrementally)

                // row upper bounds
                for (auto i = 0; i < nweights; i++)
                {
                    model.setRowUpper(i, max_weight);
                    model.setRowUpper(nweights + i, -min_weight);
                }

                // first column of the matrix
                vector<double> newcol(nweights * 2);
                for (auto j = 0; j < nweights; j++)
                {
                    newcol[j]            = fabs(EV(j, 0)) < min_ev_val ? 0.0 : EV(j, 0);
                    newcol[nweights + j] = -newcol[j];
                }
                model.addColumn(2 * nweights, &col_idx[0], &newcol[0], -COIN_DBL_MAX, COIN_DBL_MAX);

                for (auto i = 2; i <= EV.cols(); i++)        // expand from 2 eigenvectors to all, one at a time
                {
                    // add another column
                    for (auto j = 0; j < nweights; j++)
                    {
                        newcol[j]            = fabs(EV(j, i - 1)) < min_ev_val ? 0.0 : EV(j, i - 1);
                        newcol[nweights + j] = -newcol[j];
                    }
                    model.addColumn(2 * nweights, &col_idx[0], &newcol[0], -COIN_DBL_MAX, COIN_DBL_MAX);

                    // solve
                    model.setLogLevel(0);
                    model.initialSolve();               // simplex method, ~10X faster than interior point
                    //             model.primalDual();              // interior point method

                    if (!model.isProvenPrimalInfeasible() && !model.isIterationLimitReached())
                    {
                        // copy out the solution
                        VectorX<T>    solved_weights = VectorX<T>::Zero(nweights);
                        int           ncols          = model.getNumCols();
                        const double* colSol         = model.getColSolution();
                        for (auto k = 0; k < ncols; k++)
                            solved_weights += colSol[k] * EV.col(k);

                        // check if the solution was found successfully
                        if ( (solved_weights.array() > 0.0).all() )
                        {
                            weights = solved_weights;
                            weights *= (1.0 / weights.maxCoeff());  // scale to max weight = 1
                            success = true;
                            if (verbose)
                                fprintf(stderr, "curve %d: successful linear solve from %d eigenvectors\n", curve_id, i);
                        }
                    }

                    if (success)
                        break;
                }                                               // increasing number of eigenvectors

                if (!success)
                {
                    weights = VectorX<T>::Ones(nweights);
                    if (verbose)
                        fprintf(stderr, "curve %d: linear solver could not find positive weights; setting to 1\n", curve_id);
                }
            }                                                   // else need to expand eigenspace

            return success;
        }

#endif

        // computes right hand side vector of P&T eq. 9.63 and 9.67, p. 411-412 for a curve from the
        // original input domain points
        // includes multiplication by weights
        // R is column vector of n + 1 elements, each element multiple coordinates of the input points
        void RHS(
                int         cur_dim,  // current dimension
                MatrixX<T>& N,        // matrix of basis function coefficients
                MatrixX<T>& R,        // (output) residual matrix allocated by caller
                VectorX<T>& weights,  // precomputed weights for n + 1 control points on this curve
                int         po,       // index of starting parameter
                int         co)       // index of starting domain pt in current curve
        {
            int last   = R.cols() - 1;                                  // column of range value TODO: weighing only the last column does not make much sense in the split model
            MatrixX<T> Rk(N.rows(), mfa.max_dim - mfa.min_dim + 1);     // one row for each input point
            VectorX<T> denom(N.rows());                                 // rational denomoninator for param of each input point

            for (int k = 0; k < N.rows(); k++)                          // for all input points
            {
                denom(k) = (N.row(k).cwiseProduct(weights.transpose())).sum();
#ifdef UNCLAMPED_KNOTS
                if (denom(k) == 0.0)
                    denom(k) = 1.0;
#endif
                Rk.row(k) = domain.block(co + k * mfa.ds[cur_dim], mfa.min_dim, 1, mfa.max_dim - mfa.min_dim + 1);
            }

#ifdef WEIGH_ALL_DIMS                               // weigh all dimensions
            // compute the matrix R (one row for each control point)
            for (int i = 0; i < N.cols(); i++)
                for (int j = 0; j < R.cols(); j++)
                    // using array() for element-wise multiplication, which is what we want (not matrix mult.)
                    R(i, j) =
                        (N.col(i).array() *                 // ith basis functions for input pts
                         weights(i) / denom.array() *       // rationalized
                         Rk.col(j).array()).sum();          // input points
#else                                               // don't weigh domain coordinate (only range)
            // debug
//             cerr << "N in RHS:\n" << N << endl;
//             cerr << "weights in RHS:\n" << weights << endl;
//             cerr << "denom in RHS:\n" << denom << endl;

            // compute the matrix R (one row for each control point)
            for (int i = 0; i < N.cols(); i++)
            {
                // using array() for element-wise multiplication, which is what we want (not matrix mult.)
                for (int j = 0; j < R.cols() - 1; j++)
                    R(i, j) =
                        (N.col(i).array() *                 // ith basis functions for input pts
                         Rk.col(j).array()).sum();          // input points
                R(i, last) =
                    (N.col(i).array() *                     // ith basis functions for input pts
                     weights(i) / denom.array() *           // rationalized
                     Rk.col(last).array()).sum();           // input points
            }
#endif

            // debug
//             cerr << "R:\n" << R << endl;
        }


        // computes right hand side vector of P&T eq. 9.63 and 9.67, p. 411-412 for a curve from a
        // new set of input points, not the default input domain
        // includes multiplication by weights
        // R is column vector of n + 1 elements, each element multiple coordinates of the input points
        void RHS(
                int         cur_dim,  // current dimension
                MatrixX<T>& in_pts,   // input points (not the default domain stored in the mfa)
                MatrixX<T>& N,        // matrix of basis function coefficients
                MatrixX<T>& R,        // (output) residual matrix allocated by caller
                VectorX<T>& weights,  // precomputed weights for n + 1 control points on this curve
                int         po,       // index of starting parameter
                int         co,       // index of starting input pt in current curve
                int         cs)       // stride of input pts in current curve
        {
            int last   = R.cols() - 1;                                  // column of range value TODO: weighing only the last column does not make much sense in the split model
            MatrixX<T> Rk(N.rows(), mfa.max_dim - mfa.min_dim + 1);     // one row for each input point
            VectorX<T> denom(N.rows());                                 // rational denomoninator for param of each input point

            for (int k = 0; k < N.rows(); k++)
            {
                denom(k) = (N.row(k).cwiseProduct(weights.transpose())).sum();
                Rk.row(k) = in_pts.row(co + k * cs);
            }

#ifdef WEIGH_ALL_DIMS                               // weigh all dimensions
            // compute the matrix R (one row for each control point)
            for (int i = 0; i < N.cols(); i++)
                for (int j = 0; j < R.cols(); j++)
                    // using array() for element-wise multiplication, which is what we want (not matrix mult.)
                    R(i, j) =
                        (N.col(i).array() *                 // ith basis functions for input pts
                         weights(i) / denom.array() *       // rationalized
                         Rk.col(j).array()).sum();          // input points
#else                                               // don't weigh domain coordinate (only range)
            // compute the matrix R (one row for each control point)
            for (int i = 0; i < N.cols(); i++)
            {
                // using array() for element-wise multiplication, which is what we want (not matrix mult.)
                for (int j = 0; j < R.cols() - 1; j++)
                    R(i, j) =
                        (N.col(i).array() *                 // ith basis functions for input pts
                         Rk.col(j).array()).sum();          // input points
                R(i, last) =
                    (N.col(i).array() *                     // ith basis functions for input pts
                     weights(i) / denom.array() *           // rationalized
                     Rk.col(last).array()).sum();           // input points
            }
#endif

            // debug
            //     cerr << "R:\n" << R << endl;
        }

        // Checks quantities needed for approximation
        void Quants(
                VectorXi& n,          // (output) number of control point spans in each dim
                VectorXi& m)          // (output) number of input data point spans in each dim
        {
            if (mfa.p.size() != mfa.ndom_pts.size())
            {
                fprintf(stderr, "Error: Encode() size of p must equal size of ndom_pts\n");
                exit(1);
            }
            for (size_t i = 0; i < mfa.p.size(); i++)
            {
                // TODO: hard-coded for one tensor
                if (mfa.tmesh.tensor_prods[0].nctrl_pts(i) <= mfa.p(i))
                {
                    fprintf(stderr, "Error: Encode() number of control points in dimension %ld "
                            "must be at least p + 1 for dimension %ld\n", i, i);
                    exit(1);
                }
                // TODO: hard-coded for one tensor
                if (mfa.tmesh.tensor_prods[0].nctrl_pts(i) > mfa.ndom_pts(i))
                {
                    // TODO: hard-coded for one tensor
                    fprintf(stderr, "Warning: Encode() number of control points (%d) in dimension %ld "
                            "exceeds number of input data points (%d) in dimension %ld. "
                            "Technically, this is not an error, but it could be a sign something is wrong and "
                            "probably not desired if you want compression. You may not be able to get the "
                            "desired error limit and compression simultaneously. Try increasing error limit?\n",
                            mfa.tmesh.tensor_prods[0].nctrl_pts(i), i, mfa.ndom_pts(i), i);
                }
            }

            n.resize(mfa.p.size());
            m.resize(mfa.p.size());
            for (size_t i = 0; i < mfa.p.size(); i++)
            {
                // TODO: hard-coded for one tensor
                n(i)        =  mfa.tmesh.tensor_prods[0].nctrl_pts(i) - 1;
                m(i)        =  mfa.ndom_pts(i)  - 1;
            }
        }

        // solves for one curve of control points
        void CtrlCurve(
                MatrixX<T>&         N,                  // basis functions for current dimension
                MatrixX<T>&         NtN,                // Nt * N
                MatrixX<T>&         R,                  // residual matrix for current dimension and curve
                MatrixX<T>&         P,                  // solved points for current dimension and curve
                size_t              k,                  // current dimension
                size_t              co,                 // starting ofst for reading domain pts
                size_t              cs,                 // stride for reading domain points
                size_t              to,                 // starting ofst for writing control pts
                MatrixX<T>&         temp_ctrl0,         // first temporary control points buffer
                MatrixX<T>&         temp_ctrl1,         // second temporary control points buffer
                int                 curve_id,           // debugging
                TensorProduct<T>&   tensor,             // tensor product containing result
                bool                weighted = true)    // solve for and use weights
        {
            // solve for weights
            // TODO: avoid copying into Q by passing temp_ctrl0, temp_ctrl1, co, cs to Weights()
            // TODO: check that this is right, using co and cs for copying control points and domain points
            MatrixX<T> Q;
            Q.resize(mfa.ndom_pts(k), tensor.ctrl_pts.cols());
            if (k == 0)
            {
                for (auto i = 0; i < mfa.ndom_pts(k); i++)
                    Q.row(i) = domain.block(co + i * cs, mfa.min_dim, 1, mfa.max_dim - mfa.min_dim + 1);
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

            VectorX<T> weights = VectorX<T>::Ones(N.cols());

#ifndef MFA_NO_WEIGHTS

            if (weighted)
                if (k == mfa.dom_dim - 1)                      // only during last dimension of separable iteration over dimensions
                    Weights(k, Q, N, NtN, curve_id, weights);   // solve for weights

#endif

            // compute R
            // first dimension reads from domain
            // subsequent dims alternate reading temp_ctrl0 and temp_ctrl1
            // even dim reads temp_ctrl1, odd dim reads temp_ctrl0; opposite of writing order
            // because what was written in the previous dimension is read in the current one

            // TODO: move mfa.po[k] to tensor

            if (k == 0)
                RHS(k, N, R, weights, mfa.po[k], co);                 // input points = default domain
            else if (k % 2)
                RHS(k, temp_ctrl0, N, R, weights, mfa.po[k], co, cs); // input points = temp_ctrl0
            else
                RHS(k, temp_ctrl1, N, R, weights, mfa.po[k], co, cs); // input points = temp_ctrl1

            // rationalize NtN, ie, weigh the basis function coefficients
            MatrixX<T> NtN_rat = NtN;
            mfa.Rationalize(k, weights, N, NtN_rat);

            // solve for P
#ifdef WEIGH_ALL_DIMS                                   // weigh all dimensions
            P = NtN_rat.ldlt().solve(R);
#else                                                   // don't weigh domain coordinate (only range)
            // TODO: avoid 2 solves?
            MatrixX<T> P2(P.rows(), P.cols());
            P = NtN.ldlt().solve(R);                            // nonrational domain coordinates
            P2 = NtN_rat.ldlt().solve(R);                       // rational range coordinate
            for (auto i = 0; i < P.rows(); i++)
                P(i, P.cols() - 1) = P2(i, P.cols() - 1);
#endif

            // append points from P to control points that will become inputs for next dimension
            // TODO: any way to avoid this?
            CopyCtrl(P, k, co, cs, to, tensor, temp_ctrl0, temp_ctrl1);

            // copy weights of final dimension to mfa
            if (k == mfa.dom_dim - 1)
            {
                for (auto i = 0; i < weights.size(); i++)
                    tensor.weights(to + i * cs) = weights(i);
            }
        }

        // append solved control points from P to become inputs for next dimension
        // TODO: any way to avoid this copy?
        // last dimension gets copied to final control points
        // previous dimensions get copied to alternating double buffers
        void CopyCtrl(
                MatrixX<T>&         P,              // solved points for current dimension and curve
                int                 k,              // current dimension
                size_t              co,             // starting offset for reading domain points
                size_t              cs,             // stride for reading domain points
                size_t              to,             // starting offset for writing control points
                TensorProduct<T>&   tensor,         // tensor product containing result
                MatrixX<T>&         temp_ctrl0,     // first temporary control points buffer
                MatrixX<T>&         temp_ctrl1)     // second temporary control points buffer
        {
            int ndims = mfa.ndom_pts.size();    // number of domain dimensions

            // if there is only one dim, copy straight to output
            if (ndims == 1)
            {
                for (int i = 0; i < P.rows(); i++)
                    tensor.ctrl_pts.row(to + i * cs) = P.row(i);
            }
            // first dim copied from domain to temp_ctrl0
            else if (k == 0)
            {
                for (int i = 0; i < P.rows(); i++)
                    temp_ctrl0.row(to + i * cs) = P.row(i);
            }
            // even numbered dims (but not the last one) copied from P to temp_ctrl0
            else if (k % 2 == 0 && k < ndims - 1)
            {
                for (int i = 0; i < P.rows(); i++)
                    temp_ctrl0.row(to + i * cs) = P.row(i);
            }
            // odd numbered dims (but not the last one) copied from P to temp_ctrl1
            else if (k % 2 == 1 && k < ndims - 1)
            {
                for (int i = 0; i < P.rows(); i++)
                    temp_ctrl1.row(to + i * cs) = P.row(i);
            }
            // final dim if even is copied from temp_ctrl1 to ctrl_pts
            else if (k == ndims - 1 && k % 2 == 0)
            {
                for (int i = 0; i < P.rows(); i++)
                    tensor.ctrl_pts.row(to + i * cs) = P.row(i);
            }
            // final dim if odd is copied from temp_ctrl0 to ctrl_pts
            else if (k == ndims - 1 && k % 2 == 1)
            {
                for (int i = 0; i < P.rows(); i++)
                    tensor.ctrl_pts.row(to + i * cs) = P.row(i);
            }
        }

        // computes new knots to be inserted into a curve
        // for each current knot span where the error is greater than the limit, finds the domain point
        // where the error is greatest and adds the knot at that parameter value
        //
        // this version takes a set of control points as input instead of mfa.ctrl_pts
        int ErrorCurve(
                size_t         k,             // current dimension
                size_t         co,            // starting ofst for reading domain pts
                MatrixX<T>&    ctrl_pts,      // control points
                VectorX<T>&    weights,       // weights associated with control points
                VectorX<T>&    extents,       // extents in each dimension, for normalizing error (size 0 means do not normalize)
                set<int>&      err_spans,     // spans with error greater than err_limit
                T              err_limit)     // max allowable error
        {
            mfa::Decoder<T> decoder(mfa, verbose);
            int pt_dim = mfa.tmesh.tensor_prods[0].ctrl_pts.cols();     // control point dimensonality
            VectorX<T> cpt(pt_dim);                                     // decoded curve point
            int nerr = 0;                                               // number of points with error greater than err_limit
            int span = mfa.p[k];                                        // current knot span of the domain point being checked
            if (!extents.size())
                extents = VectorX<T>::Ones(domain.cols());

            for (auto i = 0; i < mfa.ndom_pts[k]; i++)      // all domain points in the curve
            {
                while (mfa.tensor.all_knots[k][span + 1] < 1.0 && mfa.tensor.all_knots[k][span + 1] <= mfa.params(mfa.po[k] + i))
                    span++;

                decoder.CurvePt(k, mfa.params(mfa.po[k] + i), ctrl_pts, weights, cpt);


                // error
                T max_err = 0.0;
                for (auto j = 0; j < mfa.max_dim - mfa.min_dim + 1; j++)
                {
                    T err = fabs(cpt(j) - domain(co + i * mfa.ds[k], mfa.min_dim + j)) / extents(mfa.min_dim + j);
                    max_err = err > max_err ? err : max_err;
                }

                if (max_err > err_limit)
                {
                    // don't duplicate spans
                    set<int>::iterator it = err_spans.find(span);
                    if (!err_spans.size() || it == err_spans.end())
                    {
                        // ensure there would be a domain point in both halves of the span if it were split
                        bool split_left = false;
                        for (auto j = i; mfa.params(mfa.po[k] + j) >= mfa.tensor.all_knots[k][span]; j--)
                            if (mfa.params(mfa.po[k] + j) < (mfa.tensor.all_knots[k][span] + mfa.tensor.all_knots[k][span + 1]) / 2.0)
                            {
                                split_left = true;
                                break;
                            }
                        bool split_right = false;
                        for (auto j = i; mfa.params(mfa.po[k] + j) < mfa.tensor.all_knots[k][span + 1]; j++)
                            if (mfa.params(mfa.po[k] + j) >= (mfa.tensor.all_knots[k][span] + mfa.tensor.all_knots[k][span + 1]) / 2.0)
                            {
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
                }
            }

            return nerr;
        }

        // encodes at full dimensionality and decodes at full dimensionality
        // decodes full-d points in each knot span and adds new knot spans where error > err_limit
        // returns true if done, ie, no knots are inserted
        bool NewKnots_full(
                VectorXi&      nnew_knots,                  // number of new knots in each dim
                vector<T>&     new_knots,                   // new knots (1st dim changes fastest)
                T              err_limit,                   // max allowable error
                VectorX<T>&    extents,                     // extents in each dimension, for normalizing error (size 0 means do not normalize)
                int            iter)                        // iteration number of caller (for debugging)
        {
            if (!extents.size())
                extents = VectorX<T>::Ones(mfa.ctrl_pts.cols());

            // resize control points and weights
            mfa.ctrl_pts.resize(mfa.tot_nctrl, mfa.max_dim - mfa.min_dim + 1);
            mfa.weights = VectorX<T>::Ones(mfa.ctrl_pts.rows());

            // full n-d encoding
            Encode();

            // find new knots
            nnew_knots = VectorXi::Zero(mfa.p.size());
            new_knots.resize(0);
            mfa::NewKnots<T> nk(mfa);
            bool done = nk.ErrorSpans(nnew_knots, new_knots, err_limit, extents, iter);

            return done;
        }

        // 1d encoding and 1d decoding
        // adds knots error spans from all curves in all directions (into a set)
        // adds knots in middles of spans that have error higher than the limit
        // returns true if done, ie, no knots are inserted
        bool NewKnots_curve(
                VectorXi&      nnew_knots,                       // number of new knots in each dim
                vector<T>&     new_knots,                        // new knots (1st dim changes fastest)
                T              err_limit,                        // max allowable error
                VectorX<T>&    extents,                          // extents in each dimension, for normalizing error (size 0 means do not normalize)
                int            iter)                             // iteration number of caller (for debugging)
        {
            // check and assign main quantities
            int  ndims  = mfa.ndom_pts.size();                          // number of domain dimensions
            int  pt_dim = mfa.tmesh.tensor_prods[0].ctrl_pts.cols();    // control point dimensonality
            // TODO: hard-coded for one tensor
            VectorXi n = mfa.tmesh.tensor_prods[0].nctrl_pts - VectorXi::Ones(ndims); // number of control point spans in each domain dim
            VectorXi m = mfa.ndom_pts  - VectorXi::Ones(ndims); // number of input data point spans in each domain dim
            nnew_knots = VectorXi::Zero(mfa.p.size());
            new_knots.resize(0);

            // resize control points and weights
            mfa.ctrl_pts.resize(mfa.tot_nctrl, mfa.max_dim - mfa.min_dim + 1);
            mfa.weights = VectorX<T>::Ones(mfa.ctrl_pts.rows());

            for (size_t k = 0; k < ndims; k++)              // for all domain dimensions
            {
                // for now set weights to 1, TODO: get weights from elsewhere
                // NB: weights are for all n + 1 control points, not just the n -1 interior ones
                VectorX<T> weights = VectorX<T>::Ones(n(k) + 1);

                // temporary control points for one curve
                // TODO: hard-coded for one tensor
                MatrixX<T> temp_ctrl = MatrixX<T>::Zero(mfa.tmesh.tensor_prods[0].nctrl_pts(k), pt_dim);

                // error spans for one curve and for worst curve
                set<int> err_spans;

                // maximum number of domain points with error greater than err_limit and their curves
                size_t max_nerr     =  0;

                // N is a matrix of (m + 1) x (n + 1) scalars that are the basis function coefficients
                //  _                          _
                // |  N_0(u[0])   ... N_n(u[0]) |
                // |     ...      ...      ...  |
                // |  N_0(u[m])   ... N_n(u[m]) |
                //  -                          -
                // TODO: N is going to be very sparse when it is large: switch to sparse representation
                // N has semibandwidth < p  nonzero entries across diagonal
                MatrixX<T> N = MatrixX<T>::Zero(m(k) + 1, n(k) + 1); // coefficients matrix

                for (int i = 0; i < N.rows(); i++)            // the rows of N
                {
                    // TODO: hard-coded for single tensor
                    int span = mfa.FindSpan(k, mfa.params(mfa.po[k] + i), mfa.tmesh.tensor_prods[0]);
                    mfa.BasisFuns(k, mfa.params(mfa.po[k] + i), span, N, i);
                }

                // TODO: NtN is going to be very sparse when it is large: switch to sparse representation
                // NtN has semibandwidth < p + 1 nonzero entries across diagonal
                MatrixX<T> NtN = N.transpose() * N;

                // R is the right hand side needed for solving NtN * P = R
                MatrixX<T> R(N.cols(), pt_dim);

                // P are the unknown interior control points and the solution to NtN * P = R
                // NtN is positive definite -> do not need pivoting
                // TODO: use a common representation for P and ctrl_pts to avoid copying
                MatrixX<T> P(N.cols(), pt_dim);

                size_t ncurves   = domain.rows() / mfa.ndom_pts(k); // number of curves in this dimension
                int nsame_steps  = 0;                                   // number of steps with same number of erroneous points
                int n_step_sizes = 0;                                   // number of step sizes so far

                // starting step size over curves
                size_t s0 = ncurves / 2 > 0 ? ncurves / 2 : 1;

                // debug, only one step size s=1
                //         s0 = 1;

                for (size_t s = s0; s >= 1 && ncurves / s < max_num_curves; s /= 2)        // for all step sizes over curves
                {
                    bool new_max_nerr = false;                      // this step size changed the max_nerr

                    for (size_t j = 0; j < ncurves; j++)            // for all the curves in this dimension
                    {
                        // each time the step changes, shift start of s-th curves by one (by subtracting
                        // n_step-sizes below)
                        if (j >= n_step_sizes && (j - n_step_sizes) % s == 0)   // this is one of the s-th curves; compute it
                        {
                            // compute R from input domain points
                            RHS(k, N, R, weights, mfa.ko[k], mfa.po[k], mfa.co[k][j]);

                            // rationalize NtN
                            MatrixX<T> NtN_rat = NtN;
                            mfa.Rationalize(k, weights, N, NtN_rat);

                            // solve for P
#ifdef WEIGH_ALL_DIMS                                           // weigh all dimensions
                            P = NtN_rat.ldlt().solve(R);
#else                                                           // don't weigh domain coordinate (only range)
                            // TODO: avoid 2 solves?
                            MatrixX<T> P2(P.rows(), P.cols());
                            P = NtN.ldlt().solve(R);                            // nonrational domain coordinates
                            P2 = NtN_rat.ldlt().solve(R);                       // rational range coordinate
                            for (auto i = 0; i < P.rows(); i++)
                                P(i, P.cols() - 1) = P2(i, P.cols() - 1);
#endif

                            // compute the error on the curve (number of input points with error > err_limit)
                            size_t nerr = ErrorCurve(k, mfa.co[k][j], P, weights, extents, err_spans, err_limit);

                            if (nerr > max_nerr)
                            {
                                max_nerr     = nerr;
                                new_max_nerr = true;
                            }
                        }
                    }                                               // curves in this dimension

                    // stop refining step if no change
                    if (max_nerr && !new_max_nerr)
                        nsame_steps++;
                    if (nsame_steps == 2)
                        break;

                    n_step_sizes++;
                }                                                   // step sizes over curves

                // free R, NtN, and P
                R.resize(0, 0);
                NtN.resize(0, 0);
                P.resize(0, 0);

                // add new knots in the middle of spans with errors
                nnew_knots(k) = err_spans.size();
                auto old_size = new_knots.size();
                new_knots.resize(old_size + err_spans.size());      // existing values are preserved
                size_t i = 0;                                       // index into new_knots
                for (set<int>::iterator it = err_spans.begin(); it != err_spans.end(); ++it)
                {
                    // debug
                    // TODO: hard-coded for one tensor
                    assert(*it < mfa.tmesh.tensor_prods[0].nctrl_pts[k]);                          // not trying to go beyond the last span

                    new_knots[old_size + i] = (mfa.tensor.all_knots[k][*it] + mfa.tensor.all_knots[k][*it + 1]) / 2.0;
                    i++;
                }

                // print progress
                //         fprintf(stderr, "\rdimension %ld of %d encoded\n", k + 1, ndims);
            }                                                      // domain dimensions

            mfa::NewKnots<T> nk(mfa);
            nk.InsertKnots(nnew_knots, new_knots);

            // debug
//             cerr << "\nnnew_knots:\n" << nnew_knots << endl;
//             cerr << "new_knots:\n"  << endl;
//             for (auto i = 0; i < new_knots.size(); i++)
//                 cerr << new_knots[i] << endl;
//             cerr << endl;

            return(nnew_knots.sum() ? 0 : 1);
        }

        template <typename>
        friend class NewKnots;

        MatrixX<T>&     domain;                         // input points
        MFA_Data<T>&    mfa;                            // the mfa object
        int             verbose;                        // output level
        size_t          max_num_curves;                 // max num. curves per dimension to check in curve version
    };
}

#endif
