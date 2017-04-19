//--------------------------------------------------------------
// decoder object
// ref: [P&T] Piegl & Tiller, The NURBS Book, 1995
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#include <mfa/mfa.hpp>
#include <mfa/decode.hpp>

#include <Eigen/Dense>

#include <vector>
#include <iostream>

using namespace std;

mfa::
Decoder::
Decoder(MFA& mfa_) :
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
    co(mfa_.co),
    knot_spans(mfa_.knot_spans),
    ndone_knot_spans(mfa_.ndone_knot_spans)
{
    // ensure that encoding was already done
    if (!p.size()         ||
        !ndom_pts.size()  ||
        !nctrl_pts.size() ||
        !domain.size()    ||
        !params.size()    ||
        !ctrl_pts.size()  ||
        !knots.size())
    {
        fprintf(stderr, "Decoder() error: Attempting to decode before encoding.\n");
        exit(0);
    }
}

// computes error in knot spans
// marks the knot spans that are done (error <= max_error in the entire span)
void
mfa::
Decoder::
ErrorSpans(float err_limit)         // max allowable error
{
    // spans that have already been split in this round (to prevent splitting twice)
    vector<bool> split_spans(knot_spans.size());    // intialized to false by default

    // debug
    fprintf(stderr, "%ld knot spans\n\n", knot_spans.size());
    fprintf(stderr, "not done spans:\n\n");

    for (auto i = 0; i < knot_spans.size(); i++)                  // knot spans
    {
        // TODO: if we swap done spans to the back, then break instead of continue
        if (knot_spans[i].done)
            continue;

        size_t nspan_pts = 1;                                   // number of domain points in the span
        for (auto k = 0; k < p.size(); k++)
            nspan_pts *= (knot_spans[i].max_param_ijk(k) - knot_spans[i].min_param_ijk(k) + 1);

        VectorXi p_ijk = knot_spans[i].min_param_ijk;           // indices of current parameter in the span
        VectorXf param(p.size());                               // value of current parameter
        bool span_done = true;                                  // span is done until error > err_limit
        // TODO:  consider binary search of the points in the span?
        // (error likely to be higher in the center of the span?)
        for (auto j = 0; j < nspan_pts; j++)                    // parameters in the span
        {
            for (auto k = 0; k < p.size(); k++)
                param(k) = params(po[k] + p_ijk(k));

            // approximate the point and measure error
            size_t idx;
            mfa.ijk2idx(p_ijk, idx);
            VectorXf cpt(ctrl_pts.cols());       // approximated point
            VolPt(param, cpt);
            float err = fabs(mfa.NormalDistance(cpt, idx)) / dom_range;     // normalized by data range

            // span is not done
            if (err > err_limit)
            {
                span_done = false;
                break;
            }


            // increment param ijk
            for (auto k = 0; k < p.size(); k++)                 // dimensions in the parameter
            {
                if (p_ijk(k) < knot_spans[i].max_param_ijk(k))
                {
                    p_ijk(k)++;
                    break;
                }
                else
                    p_ijk(k) = knot_spans[i].min_param_ijk(k);
            }                                                   // dimension in parameter
        }                                                       // parameters in the span

        // span is done
        // TODO: swap done span to the end
        if (span_done)
        {
            knot_spans[i].done = true;
            ndone_knot_spans++;
        }
    }                                                           // knot spans

    // debug
    fprintf(stderr, "number done knot spans = %ld out of %ld\n", ndone_knot_spans, knot_spans.size());
    for (auto i = 0; i < knot_spans.size(); i++)                  // knot spans
    {
        cerr <<
            "span_idx="          << i                           <<
            "\nmin_knot_ijk:\n"  << knot_spans[i].min_knot_ijk  <<
            "\nmax_knot_ijk:\n"  << knot_spans[i].max_knot_ijk  <<
            "\nmin_knot:\n"      << knot_spans[i].min_knot      <<
            "\nmax_knot:\n"      << knot_spans[i].max_knot      <<
            "\nmin_param_ijk:\n" << knot_spans[i].min_param_ijk <<
            "\nmax_param_ijk:\n" << knot_spans[i].max_param_ijk <<
            "\nmin_param:\n"     << knot_spans[i].min_param     <<
            "\nmax_param:\n"     << knot_spans[i].max_param     <<
            "\n"                 << endl;
    }

    // split spans that are not done
    auto norig_spans = knot_spans.size();
    for (auto i = 0; i < norig_spans; i++)                  // knot spans
        if (!knot_spans[i].done && !split_spans[i])
        {
            // debug
            fprintf(stderr, "\n*** calling SplitSpan(%d) ***\n\n", i);

            SplitSpan(i, split_spans);
        }

    // debug
    fprintf(stderr, "number knot spans after splitting= %ld\n", knot_spans.size());
    for (auto i = 0; i < knot_spans.size(); i++)                  // knot spans
    {
        cerr <<
            "span_idx="          << i                           <<
            "\nmin_knot_ijk:\n"  << knot_spans[i].min_knot_ijk  <<
            "\nmax_knot_ijk:\n"  << knot_spans[i].max_knot_ijk  <<
            "\nmin_knot:\n"      << knot_spans[i].min_knot      <<
            "\nmax_knot:\n"      << knot_spans[i].max_knot      <<
            "\nmin_param_ijk:\n" << knot_spans[i].min_param_ijk <<
            "\nmax_param_ijk:\n" << knot_spans[i].max_param_ijk <<
            "\nmin_param:\n"     << knot_spans[i].min_param     <<
            "\nmax_param:\n"     << knot_spans[i].max_param     <<
            "\n"                 << endl;
    }
}

// splits a knot span into two
void
mfa::
Decoder::
SplitSpan(size_t      si,               // id of span to split
        vector<bool>& split_spans)      // spans that were split already in this round, don't split these again
{
    // check if span can be split (both halves would have domain points in its range)
    // if not, check other split directions
    int sd = knot_spans[si].last_split_dim;             // new split dimension
    float new_knot;                                     // new knot value in the split dimension
    size_t k;                                           // dimension
    for (k = 0; k < p.size(); k++)
    {
        sd       = (sd + 1) % p.size();
        new_knot = (knot_spans[si].min_knot(sd) + knot_spans[si].max_knot(sd)) / 2;
        if (params(po[sd] + knot_spans[si].min_param_ijk(sd)) < new_knot &&
                params(po[sd] + knot_spans[si].max_param_ijk(sd)) > new_knot)
            break;
    }

    if (k == p.size())                                  // a split direction could not be found
    {
        knot_spans[si].done = true;
        split_spans[si]     = true;
        return;
    }


    // find all spans with the same min_knot_ijk as the span to be split and that are not done yet
    // those will be split too (NB, in the same dimension as the original span to be split)
    for (auto j = 0; j < split_spans.size(); j++)       // original number of spans in this round
    {
        if (knot_spans[j].done || split_spans[j] || knot_spans[j].min_knot_ijk(sd) != knot_spans[si].min_knot_ijk(sd))
            continue;

        // copy span to the back
        knot_spans.push_back(knot_spans[j]);

        // modify old span
        auto pi = knot_spans[j].min_param_ijk(sd);          // one coordinate of ijk index into params
        if (params(po[sd] + pi) < new_knot)                 // at least one param (domain pt) in the span
        {
            while (params(po[sd] + pi) < new_knot)          // pi - 1 = max_param_ijk(sd) in the modified span
                pi++;
            knot_spans[j].last_split_dim    = sd;
            knot_spans[j].max_knot(sd)      = new_knot;
            knot_spans[j].max_param_ijk(sd) = pi - 1;
            knot_spans[j].max_param(sd)     = params(po[sd] + pi - 1);

            // modify new span
            knot_spans.back().last_split_dim     = -1;
            knot_spans.back().min_knot(sd)       = new_knot;
            knot_spans.back().min_param_ijk(sd)  = pi;
            knot_spans.back().min_param(sd)      = params(po[sd] + pi);
            knot_spans.back().min_knot_ijk(sd)++;

            split_spans[j] = true;
        }
    }

    // increment min and max knot ijk for any knots after the inserted one
    for (auto j = 0; j < knot_spans.size(); j++)
    {
        if (knot_spans[j].min_knot(sd) > knot_spans[si].max_knot(sd))
            knot_spans[j].min_knot_ijk(sd)++;
        if (knot_spans[j].max_knot(sd) > knot_spans[si].max_knot(sd))
            knot_spans[j].max_knot_ijk(sd)++;
    }

    // debug
//     cerr << "\n\nAfter splitting spans:" << endl;
//     for (auto i = 0; i < knot_spans.size(); i++)
//     {
//         if (knot_spans[i].done)
//             continue;
//         cerr <<
//             "span_idx="          << i                            <<
//             "\nmin_knot_ijk:\n"  << knot_spans[i].min_knot_ijk  <<
//             "\nmax_knot_ijk:\n"  << knot_spans[i].max_knot_ijk  <<
//             "\nmin_knot:\n"      << knot_spans[i].min_knot      <<
//             "\nmax_knot:\n"      << knot_spans[i].max_knot      <<
//             "\nmin_param_ijk:\n" << knot_spans[i].min_param_ijk  <<
//             "\nmax_param_ijk:\n" << knot_spans[i].max_param_ijk  <<
//             "\nmin_param:\n"     << knot_spans[i].min_param      <<
//             "\nmax_param:\n"     << knot_spans[i].max_param      <<
//             "\n"                 << endl;
//     }
}

// computes approximated points from a given set of domain points and an n-d NURBS volume
// P&T eq. 9.77, p. 424
// this version recomputes parameter values of input points and
// recomputes basis functions rather than taking them as an input
// this version also assumes weights = 1; no division by weight is done
// assumes all vectors have been correctly resized by the caller
void
mfa::
Decoder::
Decode(MatrixXf& approx)                 // pts in approximated volume (1st dim. changes fastest)
{
    vector<size_t> iter(p.size(), 0);    // parameter index (iteration count) in current dim.
    vector<size_t> ofst(p.size(), 0);    // start of current dim in linearized params

    for (size_t i = 0; i < p.size() - 1; i++)
        ofst[i + 1] = ofst[i] + ndom_pts(i);

    // eigen frees following temp vectors when leaving scope
    VectorXf dpt(domain.cols());         // original data point
    VectorXf cpt(ctrl_pts.cols());       // approximated point
    VectorXf d(domain.cols());           // apt - dpt
    VectorXf param(p.size());            // parameters for one point
    for (size_t i = 0; i < domain.rows(); i++)
    {
        // debug
        // cerr << "input point:\n" << domain.row(i) << endl;

        // extract parameter vector for one input point from the linearized vector of all params
        for (size_t j = 0; j < p.size(); j++)
            param(j) = params(iter[j] + ofst[j]);

        // compute approximated point for this parameter vector
        VolPt(param, cpt);

        // update the indices in the linearized vector of all params for next input point
        for (size_t j = 0; j < p.size(); j++)
        {
            if (iter[j] < ndom_pts(j) - 1)
            {
                iter[j]++;
                break;
            }
            else
                iter[j] = 0;
        }

        approx.row(i) = cpt;

        // print progress
        if (i > 0 && domain.rows() >= 100 && i % (domain.rows() / 100) == 0)
            fprintf(stderr, "\r%.0f %% decoded", (float)i / (float)(domain.rows()) * 100);
    }
    fprintf(stderr, "\r100 %% decoded\n");
}

// compute a point from a NURBS curve at a given parameter value
// algorithm 4.1, Piegl & Tiller (P&T) p.124
// this version recomputes basis functions rather than taking them as an input
// this version also assumes weights = 1; no division by weight is done
void
mfa::
Decoder::
CurvePt(int       cur_dim,                     // current dimension
        float     param,                       // parameter value of desired point
        VectorXf& out_pt)                      // (output) point
{
    int n      = (int)ctrl_pts.rows() - 1;     // number of control point spans
    int span   = mfa.FindSpan(cur_dim, param);
    MatrixXf N = MatrixXf::Zero(1, n + 1);     // basis coefficients
    mfa.BasisFuns(cur_dim, param, span, N, 0, n, 0);
    out_pt = VectorXf::Zero(ctrl_pts.cols());  // initializes and resizes

    for (int j = 0; j <= p(cur_dim); j++)
        out_pt += N(0, j + span - p(cur_dim)) * ctrl_pts.row(span - p(cur_dim) + j);

    // debug
    // cerr << "n " << n << " param " << param << " span " << span << " out_pt " << out_pt << endl;
    // cerr << " N " << N << endl;
}

// decode one curve in one dimension (e.g., y direction) from a set of previously computed
// control curves in the prior dimension (e.g., x direction)
// new curve in current dimension (e.g., y) is normal to curves in prior dimension (e.g., x) and
// is located in the prior dimension (e.g., x) at parameter value prev_param
//
// assumes the caller resized the output out_pts to the correct size, which is number of original
// domain points in the current dimension (e.g., y)
//
// currently not used but can be useful for getting a cross-section curve from a surface
// would need to be expanded to get a curve from a higher dimensional space, currently gets 1D curve
// from 2D surface
void
mfa::
Decoder::
DecodeCurve(size_t    cur_dim,    // current dimension
            float     pre_param,  // parameter value in prior dimension of the pts in the curve
            size_t    ko,         // starting offset for knots in current dim
            size_t    cur_cs,     // stride for control points in current dim
            size_t    pre_cs,     // stride for control points in prior dim
            MatrixXf& out_pts)    // output approximated pts for the curve
{
    if (cur_dim == 0 || cur_dim >= p.size())
    {
        fprintf(stderr, "Error in DecodeCurve(): "
                "cur_dim out of range (must be 1 <= cur_dim <= %ld\n",
                p.size() - 1);
        exit(0);
    }

    // get one set of knots in the prior dimension
    size_t nknots = nctrl_pts(cur_dim - 1) + p(cur_dim - 1) + 1;
    VectorXf pre_knots(nknots);
    for (size_t i = 0; i < nknots; i++)
        pre_knots(i) = knots(ko + i);

    // debug
    // cerr << "pre_knots:\n" << pre_knots << endl;

    // control points for one curve in prior dimension
    MatrixXf pre_ctrl_pts(nctrl_pts(cur_dim - 1), domain.cols());
    size_t n = 0;                            // index of current control point

    // for the desired approximated points in the current dim, same number as domain points
    for (size_t j = 0; j < ndom_pts(cur_dim); j++)
    {
        // get one curve of control points in the prior dim
        for (size_t i = 0; i < nctrl_pts(cur_dim - 1); i++)
        {
            pre_ctrl_pts.row(i) = ctrl_pts.row(n);
            n += pre_cs;
        }
        n += (cur_cs - nctrl_pts(cur_dim - 1));

        // debug
        // cerr << "pre_ctrl_pts:\n" << pre_ctrl_pts << endl;

        // get the approximated point
        VectorXf out_pt(domain.cols());
        CurvePt(cur_dim - 1, pre_param, out_pt);
        out_pts.row(j) = out_pt;

        // debug
        // cerr << "out_pt: " << out_pt << endl;
    }

    // debug
    // cerr << "out_pts:\n " << out_pts << endl;
}

// compute a point from a NURBS n-d volume at a given parameter value
// algorithm 4.3, Piegl & Tiller (P&T) p.134
// this version recomputes basis functions rather than taking them as an input
// this version also assumes weights = 1; no division by weight is done
//
// There are two types of dimensionality:
// 1. The dimensionality of the NURBS tensor product (p.size())
// (1D = NURBS curve, 2D = surface, 3D = volumem 4D = hypervolume, etc.)
// 2. The dimensionality of individual control points (ctrl_pts.cols())
// p.size() should be < ctrl_pts.cols()
void
mfa::
Decoder::
VolPt(VectorXf& param,                       // parameter value in each dim. of desired point
      VectorXf& out_pt)                      // (output) point
{
    // debug
    // cerr << "\n\nparam:\n" << param << endl;

    // check dimensionality for sanity
    assert(p.size() < ctrl_pts.cols());

    out_pt = VectorXf::Zero(ctrl_pts.cols());// initializes and resizes
    vector <MatrixXf> N(p.size());           // basis functions in each dim.
    vector<VectorXf>  temp(p.size());        // temporary point in each dim.
    vector<int>       span(p.size());        // span in each dim.
    vector<int>       n(p.size());           // number of control point spans in each dim
    vector<int>       iter(p.size());        // iteration number in each dim.
    int               tot_iters = 1;         // tot. num. iterations in flattened n-d nested loops
    vector<size_t>    ko(p.size(), 0);       // starting offset for knots in current dim
    vector<size_t>    co(p.size());          // starting offset for control points in each dim
    vector<size_t>    ct(p.size());          // relative coordinates of ctrl pt of current iteration
    vector<size_t>    cs(p.size());          // stride for next co in each dim
    VectorXf          ctrl_pt(ctrl_pts.cols()); // one control point

    for (size_t i = 0; i < p.size(); i++)    // for all dims
    {
        temp[i]    = VectorXf::Zero(ctrl_pts.cols());
        iter[i]    = 0;
        tot_iters  *= (p(i) + 1);
        n[i]       = (int)nctrl_pts(i) - 1;
        span[i]    = mfa.FindSpan(i, param(i), ko[i]);
        N[i]       = MatrixXf::Zero(1, n[i] + 1);
        mfa.BasisFuns(i, param(i), span[i], N[i], 0, n[i], 0, ko[i]);
        if (i == 0)
            cs[i] = 1;
        else
            cs[i] = cs[i - 1] * nctrl_pts(i);
        co[i] = span[i] - p(i) - ko[i];
        if (i < p.size() - 1)
            ko[i + 1] = ko[i] + n[i] + p(i) + 2; // n[i]+p(i)+2 =  number of knots in current dim.

        // debug
        // cerr << "N[" << i << "] = " << N[i] << endl;
        // fprintf(stderr, "i=%d co[i]=%d span[i]=%d p(i)=%d ko[i]=%d cs[i]=%d\n",
        //         i, co[i], span[i], p(i), ko[i], cs[i]);
    }

    for (int i = 0; i < tot_iters; i++)      // 1-d flattening all n-d nested loop computations
    {
        // compute coordinates of first control point of curve corresponding to this iteration
        // these are relative to start of the box of control points located at co
        int div = tot_iters;
        int i_temp = i;
        for (int j = p.size() - 1; j >= 0; j--)
        {
            div    /= (p(j) + 1);
            ct[j]  = i_temp / div;
            i_temp -= (ct[j] * div);
        }

        // control point linear order index
        int ctrl_idx = 0;
        int fac = 1;
        for (int j = 0; j < p.size(); j++)
        {
            ctrl_idx += (co[j] + ct[j]) * fac;
            fac *= nctrl_pts(j);
        }

        // always compute the point in the first dimension
        ctrl_pt = ctrl_pts.row(ctrl_idx);
        temp[0] += (N[0])(0, iter[0] + span[0] - p(0)) * ctrl_pt;

        // debug
//         fprintf(stderr, "1: temp[0] += N[0, %d] * ctrl_pt\n", iter[0] + span[0] - p(0));
//         cerr << "span: " << span[0] << " ctrl_pt=\n" << ctrl_pt << endl;
        // cerr << "1: i=" << i << " iter[0]=" << iter[0] << " temp[0]=\n" << temp[0] << endl;
        iter[0]++;

        // for all dimensions except last, check if span is finished
        for (size_t k = 0; k < p.size() - 1; k++)
        {
            if (iter[k] - 1 == p(k))
            {
                // compute point in next higher dimension
                temp[k + 1] += (N[k + 1])(0, iter[k + 1] + span[k + 1] - ko[k + 1] - p(k + 1)) *
                    temp[k];

                // debug
                // fprintf(stderr, "2: temp[%d] += N[%d, %d] * temp[%d]\n",
                //         k + 1, k + 1, iter[k + 1] + span[k + 1] - ko[k + 1] - p(k + 1), k);
                // fprintf(stderr, "3: temp[%d] = 0\n", k);
                // cerr << "2: i=" << i << " k=" << k << " iter[k+1]=" << iter[k+1] <<
                //     " temp[k+1]=\n" << temp[k+1] << endl;

                // reset the computation for the current dimension
                temp[k]    = VectorXf::Zero(ctrl_pts.cols());
                iter[k]    = 0;
                iter[k + 1]++;
            }
        }
    }

    out_pt = temp[p.size() - 1];

    // debug
    // fprintf(stderr, "out_pt = temp[%d]\n", p.size() - 1);
    // cerr << "out_pt:\n" << out_pt << endl;
}
