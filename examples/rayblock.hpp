//--------------------------------------------------------------
// One diy block that can handle ray integrals
//
// David Lenz
// Argonne National Laboratory
// dlenz@anl.gov
//--------------------------------------------------------------
#ifndef _MFA_RAY_BLOCK
#define _MFA_RAY_BLOCK

#include "block.hpp"

using namespace std;

template <typename T>
struct RayBlock : public Block<T>
{
    using Base = typename Block<T>::BlockBase<T>;
    using Base::dom_dim;
    using Base::pt_dim;
    using Base::core_mins;
    using Base::core_maxs;
    using Base::bounds_mins;
    using Base::bounds_maxs;
    using Base::overlaps;
    using Base::input;
    using Base::approx;
    using Base::errs;
    using Base::mfa;

    // using Block<T>::create;
    // using Block<T>::destroy;
    // using Block<T>::add;
    // using Block<T>::save;
    // using Block<T>::load;

    static
        void* create()              { return mfa::create<RayBlock>(); }

    static
        void destroy(void* b)       { mfa::destroy<RayBlock>(b); }

    static
        void add(                                   // add the block to the decomposition
            int                 gid,                // block global id
            const Bounds<T>&    core,               // block bounds without any ghost added
            const Bounds<T>&    bounds,             // block bounds including any ghost region added
            const Bounds<T>&    domain,             // global data bounds
            const RCLink<T>&    link,               // neighborhood
            diy::Master&        master,             // diy master
            int                 dom_dim,            // domain dimensionality
            int                 pt_dim,             // point dimensionality
            T                   ghost_factor = 0.0) // amount of ghost zone overlap as a factor of block size (0.0 - 1.0)
    {
        mfa::add<RayBlock, T>(gid, core, bounds, domain, link, master, dom_dim, pt_dim, ghost_factor);
    }

    static
        void save(const void* b_, diy::BinaryBuffer& bb)    { mfa::save<RayBlock, T>(b_, bb); }
    static
        void load(void* b_, diy::BinaryBuffer& bb)          { mfa::load<RayBlock, T>(b_, bb); }

    RayBlock() :
        ray_mfa(nullptr),
        ray_input(nullptr),
        ray_approx(nullptr),
        ray_errs(nullptr)
    { }

    ~RayBlock()
    {
        delete ray_mfa;
        delete ray_input;
        delete ray_approx;
        delete ray_errs;
    }

    int                 ray_dom_dim;                 // dom_dim of the extended model (i.e. dom_dim+1)
    mfa::MFA<T>         *ray_mfa;
    mfa::PointSet<T>    *ray_input;                 // input data
    mfa::PointSet<T>    *ray_approx;                // output data
    mfa::PointSet<T>    *ray_errs;                  // error field

    VectorX<T>  box_mins;
    VectorX<T>  box_maxs;   

    VectorX<T>          ray_bounds_mins;            // local domain minimum corner
    VectorX<T>          ray_bounds_maxs;            // local domain maximum corner
    VectorX<T>          ray_core_mins;              // local domain minimum corner w/o ghost
    VectorX<T>          ray_core_maxs;              // local domain maximum corner w/o ghost

    vector<T>           ray_max_errs;               // maximum (abs value) distance from input points to curve
    vector<T>           ray_sum_sq_errs;            // sum of squared errors

    void get_box_intersections(
        T alpha,
        T rho,
        T& x0,
        T& y0,
        T& x1,
        T& y1,
        const VectorX<T>& mins,
        const VectorX<T>& maxs) const
    {
        T xl = mins(0);
        T xh = maxs(0);
        T yl = mins(1);
        T yh = maxs(1);

        T yh_int = (rho - yh * sin(alpha)) / cos(alpha);
        T yl_int = (rho - yl * sin(alpha)) / cos(alpha);
        T xh_int = (rho - xh * cos(alpha)) / sin(alpha);
        T xl_int = (rho - xl * cos(alpha)) / sin(alpha);
        // T x0, x1, y0, y1;
 
        // cerr << "ia=" << ia << ", ir=" << ir << endl;
        // cerr << "rho=" << rho << ", alpha=" << alpha << endl;
        // cerr << xl_int << " " << xh_int << " " << yl_int << " " << yh_int << endl;

        // "box intersection" setup
        // start/end coordinates of the ray formed by intersecting 
        // the line with bounding box of the data
        if (alpha == 0)    // vertical lines (top to bottom)
        {
            x0 = rho;
            y0 = yh;
            x1 = rho;
            y1 = yl;
        }
        else if (sin(alpha) == 0 && alpha > 0) // vertical lines (bottom to top)
        {
            x0 = rho;
            y0 = yl;
            x1 = rho;
            y1 = yh;
        }
        else if (cos(alpha)==0) // horizontal lines
        {
            x0 = xl;
            y0 = rho;
            x1 = xh;
            y1 = rho;
        }
        else if (xl_int >= yl && xl_int <= yh)  // enter left
        {
            x0 = xl;
            y0 = xl_int;

            if (yl_int >= xl && yl_int <= xh)   // enter left, exit bottom
            {
                y1 = yl;
                x1 = yl_int;
            }
            else if (yh_int >= xl && yh_int <= xh)  // enter left, exit top
            {
                y1 = yh;
                x1 = yh_int;
            }
            else if (xh_int >= yl && xh_int <= yh)  // enter left, exit right
            {
                x1 = xh;
                y1 = xh_int;
            }
            else
            {
                cerr << "ERROR: invalid state 1" << endl;
                // cerr << "ia = " << ia << ", ir = " << ir << endl;
                exit(1);
            }
        }
        else if (yl_int >= xl && yl_int <= xh)  // enter or exit bottom
        {
            if (yh_int >= xl && yh_int <= xh)   // enter/exit top & bottom
            {
                if (sin(alpha) == 0)    // vertical line case (should have been handled above)
                {
                    cerr << "ERROR: invalid state 6" << endl;
                    x0 = yl_int;
                    y0 = yl;
                    x1 = yh_int;
                    y1 = yh;
                }
                else if (sin(alpha) == 0 && alpha > 0)     // opposite vertical line case (should have been handled above)
                {
                    cerr << "ERROR: invalid state 7" << endl;
                    x0 = yh_int;
                    y0 = yh;
                    x1 = yl_int;
                    y1 = yl;
                }
                // else if (yl_int < yh_int)   // enter bottom, exit top
                else if (alpha > 3.14159265358979/2)
                { 
                    x0 = yl_int;
                    y0 = yl;
                    x1 = yh_int;
                    y1 = yh;
                }
                // else if (yl_int > yh_int)   // enter top, exit bottom
                else if (alpha < 3.14159265358979/2)
                {
                    x0 = yh_int;
                    y0 = yh;
                    x1 = yl_int;
                    y1 = yl;
                }
                else
                {
                    cerr << "ERROR: invalid state 2" << endl;
                    // cerr << "ia = " << ia << ", ir = " << ir << endl;
                    exit(1);
                }
            }
            else if (xh_int >= yl && xh_int <= yh)  // enter bottom, exit right
            {
                x0 = yl_int;
                y0 = yl;
                x1 = xh;
                y1 = xh_int;
            }
            else
            {
                cerr << "ERROR: invalid state 3" << endl;
                // cerr << "ia = " << ia << ", ir = " << ir << endl;
                exit(1);
            }
        }
        else if (yh_int >= xl && yh_int <= xh)  // enter top (cannot be exit top b/c of cases handled previously)
        {
            if (xh_int >= yl && xh_int <= yh)   // enter top, exit right
            {
                x0 = yh_int;
                y0 = yh;
                x1 = xh;
                y1 = xh_int;
            }
            else
            {
                cerr << "ERROR: invalid state 4" << endl;
                // cerr << "ia = " << ia << ", ir = " << ir << endl;
                exit(1);
            }
        }
        else
        {
            x0 = 0;
            y0 = 0;
            x1 = 0;
            y1 = 0;
            // cerr << "ERROR: invalid state 5" << endl;
            // cerr << "ia = " << ia << ", ir = " << ir << endl;
            // exit(1);
        }
    }

    // ONLY 2d AT THE MOMENT
    // precondition: Block already contains a fully encoded MFA
    void create_ray_model(
        const       diy::Master::ProxyWithLink& cp,
        MFAInfo& mfa_info,
        DomainArgs& args,
        bool fixed_length,
        int n_samples,
        int n_rho,
        int n_alpha,
        int v_samples,
        int v_rho,
        int v_alpha)
    {
        const double pi = 3.14159265358979;
        if (n_samples == 0 || n_rho == 0 || n_alpha == 0)
        {
            cerr << "ERROR: Did not set n_samples, n_rho, or n_alpha before creating a ray model. See command line help" << endl;
            exit(1);
        }
        if (v_samples == 0 || v_rho == 0 || v_alpha == 0)
        {
            cerr << "ERROR: Did not set v_samples, v_rho, or v_alpha before creating a ray model. See command line help" << endl;
            exit(1);
        }

        ray_dom_dim = dom_dim + 1;   // new dom_dim
        int new_pd = pt_dim + 1;    // new pt_dim

        // Create a model_dims vector for the auxilliary model, and increase the dimension of geometry
        VectorXi new_mdims = mfa->model_dims();
        new_mdims[0] += 1;  

        VectorXi ndom_pts{{n_samples, n_rho, n_alpha}};
        int npts = ndom_pts.prod();

        ray_input = new mfa::PointSet<T>(ray_dom_dim, new_mdims, npts, ndom_pts);

        // extents of domain in physical space
        VectorX<T> param(dom_dim);
        VectorX<T> outpt(pt_dim);
        const T xl = bounds_mins(0);
        const T xh = bounds_maxs(0);
        const T yl = bounds_mins(1);
        const T yh = bounds_maxs(1);

        box_mins = bounds_mins.head(dom_dim);
        box_maxs = bounds_maxs.head(dom_dim);

        // TODO: make this generic
        double r_lim = 0;
        if (fixed_length)
        {
            double max_radius = max(max(abs(xl),abs(xh)), max(abs(yl),abs(yh)));
            r_lim = max_radius * 1.5;
        } 
        else
        {
            r_lim = 0.99 * xh; // HACK this only works for square domains centered at origin, and for the "box intersection" setup
        }
        double dr = r_lim * 2 / (n_rho-1);
        double da = pi / (n_alpha-1); // d_alpha; amount to rotate on each slice

        // fill ray data set
        double alpha    = 0;   // angle of rotation
        double rho      = -r_lim;
        for (int ia = 0; ia < n_alpha; ia++)
        {
            alpha = ia * da;

            for (int ir = 0; ir < n_rho; ir++)
            {
                rho = -r_lim + ir * dr;

                T x0, y0, x1, y1, span_x, span_y;
                if (fixed_length)
                {
                    // "parallel-plate setup"
                    // start/end coordinates of the ray (alpha, rho)
                    // In this setup the length of every segment (x0,y0)--(x1,y1) is constant
                    span_x = 2 * r_lim * sin(alpha);
                    span_y = 2 * r_lim * cos(alpha);
                    x0 = rho * cos(alpha) - r_lim * sin(alpha);
                    x1 = rho * cos(alpha) + r_lim * sin(alpha);
                    y0 = rho * sin(alpha) + r_lim * cos(alpha);
                    y1 = rho * sin(alpha) - r_lim * cos(alpha);
                }
                else
                {
                    get_box_intersections(alpha, rho, x0, y0, x1, y1, this->box_mins, this->box_maxs);
                    span_x = x1 - x0;
                    span_y = y1 - y0;
                }

                T dx = span_x / (n_samples-1);
                T dy = span_y / (n_samples-1);

                for (int is = 0; is < n_samples; is++)
                {
                    int idx = ia*n_rho*n_samples + ir*n_samples + is;
                    ray_input->domain(idx, 0) = (double)is / (n_samples-1);
                    ray_input->domain(idx, 1) = rho;
                    ray_input->domain(idx, 2) = alpha;

                    T x = x0 + is * dx;
                    T y = 0;

                    if (fixed_length)
                        y = y0 - is * dy;
                    else
                        y = y0 + is * dy;

                    // If this point is not in the original domain
                    if (x < xl - 1e-8 || x > xh + 1e-8 || y < yl - 1e-8 || y > yh + 1e-8)
                    {
                        if (fixed_length)  // do nothing in fixed_length setting
                        {
                            ray_input->domain(idx, ray_dom_dim) = 1000;
                            continue;
                        }
                        else                // else complain and zero-pad (this should not happen)
                        {
                            cerr << "NOT IN DOMAIN" << endl;
                            cerr << "  " << x << "\t" << y << endl;
                            ray_input->domain(idx,3) = 0;
                        }
                    }
                    else    // point is in domain, decode value from existing MFA
                    {
                        param(0) = (x - xl) / (xh - xl);
                        param(1) = (y - yl) / (yh - yl);

                        // Truncate to [0,1] in the presence of small round-off errors
                        param(0) = param(0) < 0 ? 0 : param(0);
                        param(1) = param(1) < 0 ? 0 : param(1);
                        param(0) = param(0) > 1 ? 1 : param(0);
                        param(1) = param(1) > 1 ? 1 : param(1);

                        outpt.resize(pt_dim);
                        this->mfa->Decode(param, outpt);
                        ray_input->domain.block(idx, ray_dom_dim, 1, pt_dim - dom_dim) = outpt.tail(pt_dim - dom_dim).transpose();
                    }
                }
            }
        }
        
        // Set parameters for new input
        VectorX<T> input_mins(ray_dom_dim), input_maxs(ray_dom_dim);
        if (fixed_length)
        {
            input_mins(0) = 0; input_maxs(0) = 1;
            input_mins(1) = -r_lim; input_maxs(1) = r_lim;
            input_mins(2) = 0; input_maxs(2) = pi;
        }
        cerr << "input_mins: " << input_mins << endl;
        cerr << "input_maxs: " << input_maxs << endl;

        if (fixed_length)
            ray_input->set_bounds(input_mins, input_maxs);

        ray_input->set_domain_params();

        // ------------ Creation of new MFA ------------- //
        //
        // Create a new top-level MFA
        int verbose = mfa_info.verbose && cp.master()->communicator().rank() == 0; 
        ray_mfa = new mfa::MFA<T>(ray_dom_dim, verbose);

        // Set up new geometry
        ray_mfa->AddGeometry(ray_dom_dim);

        // Set nctrl_pts, degree for variables
        VectorXi nctrl_pts(ray_dom_dim);
        VectorXi p(ray_dom_dim);
        for (auto i = 0; i< this->mfa->nvars(); i++)
        {
            int min_p = 20;
            int max_nctrl_pts = 0;
            for (int j = 0; j < dom_dim; j++)   // Loop over dimensions of original MFA
            {
                if (this->mfa->var(i).p(j) < min_p)
                    min_p = this->mfa->var(i).p(j);

                if (this->mfa->var(i).tmesh.tensor_prods[0].nctrl_pts(j) > max_nctrl_pts)
                    max_nctrl_pts = this->mfa->var(i).tmesh.tensor_prods[0].nctrl_pts(j);
            }

            // set ray model degree to minimum degree of original model
            for (auto j = 0; j < ray_dom_dim; j++)
            {
                p(j)            = min_p;
            }

            nctrl_pts(0) = v_samples;
            nctrl_pts(1) = v_rho;
            nctrl_pts(2) = v_alpha;
            // p(0) = 2;
            // p(1) = 2;
            // p(2) = 2;

            ray_mfa->AddVariable(p, nctrl_pts, 1);
        }

        // Encode ray model. 
        ray_mfa->FixedEncodeGeom(*ray_input, 0, false);
        ray_mfa->RayEncode(0, *ray_input);

        // ----------- Replace old block members with new ---------- //
        // reset block members as needed
        // dom_dim = new_dd;
        // pt_dim = new_pd;

        // // replace original mfa with ray-mfa
        // delete this->mfa;
        // this->mfa = ray_mfa;
        // ray_mfa = nullptr;

        // // replace original input with ray-model input
        // delete input;
        // input = ray_input;
        // ray_input = nullptr;

        // if (new_pd != this->mfa->pt_dim)    // sanity check
        // {
        //     cerr << "ERROR: pt_dim does not match in create_ray_model()" << endl;
        //     exit(1);
        // }

        // VectorX<T> old_bounds_mins = bounds_mins;
        // VectorX<T> old_bounds_maxs = bounds_maxs;

        // bounds_mins.resize(pt_dim);
        // bounds_maxs.resize(pt_dim);
        // bounds_mins(0) = 0;
        // bounds_maxs(0) = 1;
        // bounds_mins(1) = -r_lim;
        // bounds_maxs(1) = r_lim;
        // bounds_mins(2) = 0;
        // bounds_maxs(2) = pi;
        // for (int i = 0; i < this->mfa->pt_dim - this->mfa->geom_dim(); i++)
        // {
        //     bounds_mins(3+i) = old_bounds_mins(2+i);
        //     bounds_maxs(3+i) = old_bounds_maxs(2+i);
        // }
        // core_mins = bounds_mins.head(dom_dim);
        // core_maxs = bounds_maxs.head(dom_dim);

        ray_bounds_mins.resize(pt_dim + 1);
        ray_bounds_maxs.resize(pt_dim + 1);
        ray_bounds_mins(0) = 0;
        ray_bounds_maxs(0) = 1;
        ray_bounds_mins(1) = -r_lim;
        ray_bounds_maxs(1) = r_lim;
        ray_bounds_mins(2) = 0;
        ray_bounds_maxs(2) = pi;
        for (int i = dom_dim; i < this->mfa->pt_dim; i++)
        {
            ray_bounds_mins(i+1) = bounds_mins(i);
            ray_bounds_maxs(i+1) = bounds_maxs(i);
        }
        ray_core_mins = ray_bounds_mins.head(dom_dim+1);
        ray_core_maxs = ray_bounds_maxs.head(dom_dim+1);

        ray_max_errs.resize(ray_mfa->nvars());
        ray_sum_sq_errs.resize(ray_mfa->nvars());

        // --------- Decode (for visualization) --------- //
        // this->decode_block(cp, 0);
        // this->range_error(cp, true, true);
        cout << "Decoding Ray Model to uniform grid..." << endl;
        vector<int> grid_size{{n_samples, n_rho, n_alpha}};
        VectorXi gridpoints(3);
        gridpoints(0) = grid_size[0];
        gridpoints(1) = grid_size[1];
        gridpoints(2) = grid_size[2];
        decode_ray_block(cp);

cerr << "  ===========" << endl;
cerr << "  f(x) = sin(x) hardcoded in create_ray_model() for error computation" << endl;
cerr << "  ===========" << endl;
        // delete this->errs;
        ray_errs = new mfa::PointSet<T>(ray_dom_dim, ray_mfa->model_dims(), gridpoints.prod(), gridpoints);
        outpt = VectorX<T>::Zero(1);
        param.resize(ray_dom_dim); // n.b. this is now the NEW dom_dim
        for (int k = 0; k < grid_size[2]; k++)
        {
            for (int j = 0; j < grid_size[1]; j++)
            {
                T rh_param = (T) j / (grid_size[1]-1);
                T al_param = (T) k / (grid_size[2]-1);
                T rh = ray_input->mins(1) + (ray_input->maxs(1) - ray_input->mins(1)) * rh_param;
                T al = ray_input->mins(2) + (ray_input->maxs(2) - ray_input->mins(2)) * al_param;

                T x0, y0, x1, y1, span_x, span_y;

                // "parallel-plate setup"
                // start/end coordinates of the ray (alpha, rho)
                // In this setup the length of every segment (x0,y0)--(x1,y1) is constant
                span_x = 2 * r_lim * sin(al);
                span_y = 2 * r_lim * cos(al);
                x0 = rh * cos(al) - r_lim * sin(al);
                x1 = rh * cos(al) + r_lim * sin(al);
                y0 = rh * sin(al) + r_lim * cos(al);
                y1 = rh * sin(al) - r_lim * cos(al);

                T dx = span_x / (grid_size[0]-1);
                T dy = span_y / (grid_size[0]-1);

                for (int i = 0; i < grid_size[0]; i++)
                {
                    T t_param = (T)i / (grid_size[0]-1);
                    int idx = k*grid_size[0]*grid_size[1] + j*grid_size[0] + i;
                    T a = ray_input->mins(0) + (ray_input->maxs(0) - ray_input->mins(0)) / (grid_size[0]-1) * i;
                    
                    T x = x0 + i * dx;
                    T y = 0;

                    if (fixed_length)
                        y = y0 - i * dy;
                    else
                        y = y0 + i * dy;
                    
                    param(0) = t_param;
                    param(1) = rh_param;
                    param(2) = al_param;

                    // Truncate to [0,1] in the presence of small round-off errors
                    param(0) = param(0) < 0 ? 0 : param(0);
                    param(1) = param(1) < 0 ? 0 : param(1);
                    param(2) = param(2) < 0 ? 0 : param(2);
                    param(0) = param(0) > 1 ? 1 : param(0);
                    param(1) = param(1) > 1 ? 1 : param(1);
                    param(2) = param(2) > 1 ? 1 : param(2);

                    ray_mfa->DecodeVar(0, param, outpt);

                    T trueval = sin(x) * sin(y);

                    ray_errs->domain(idx, 0) = t_param;
                    ray_errs->domain(idx, 1) = rh;
                    ray_errs->domain(idx, 2) = al;
                    ray_errs->domain(idx, 3) = abs(trueval - outpt(0));

                    // ignore "errors" when querying outside the domain
                    if (x < xl || x > xh || y < yl || y > yh)
                        ray_errs->domain(idx, 3) = 0;
                }
            }
        }
        cout << "done." << endl;

        // Compute error metrics
        MatrixX<T>& errpts = ray_errs->domain;

        for (auto j = ray_dom_dim; j < ray_errs->pt_dim; j++)
            ray_sum_sq_errs[j - ray_dom_dim] = 0.0;
        for (auto i = 0; i < ray_errs->npts; i++)
        {
            for (auto j = ray_dom_dim; j < ray_errs->pt_dim; j++)
            {
                ray_sum_sq_errs[j - ray_dom_dim] += (errpts(i, j) * errpts(i, j));
                if ((i == 0 && j == ray_dom_dim) || errpts(i, j) > ray_max_errs[j - ray_dom_dim])
                    ray_max_errs[j - ray_dom_dim] = errpts(i, j);
            }
        }
    }

    // decode entire rotation space  at the same parameter locations as 'ray_input'
    void decode_ray_block(const diy::Master::ProxyWithLink& cp)
    {
        if (ray_approx)
        {
            cerr << "WARNING: Overwriting \"ray_approx\" pointset in RayBlock::decode_ray_block" << endl;
            delete ray_approx;
        }
        ray_approx = new mfa::PointSet<T>(ray_input->params, ray_input->model_dims());  // Set decode params from ray_input params

        ray_mfa->Decode(*ray_approx, false);
    }

    pair<T,T> dualCoords(const VectorX<T>& a, const VectorX<T>& b) const
    {
        const double pi = 3.14159265358979;

        T a_x = a(0);
        T a_y = a(1);
        T b_x = b(0);
        T b_y = b(1);

         // distance in x and y between the endpoints of the segment
        T delta_x = b_x - a_x;
        T delta_y = b_y - a_y;


        T alpha = -1;
        T rho = 0;

        if (a_x == b_x)
        {
            alpha = 0;
            rho = a_x;
        }
        else
        {
            T m = (b_y-a_y)/(b_x-a_x);
            alpha = pi/2 - atan(-m);            // acot(x) = pi/2 - atan(x)
            rho = (a_y - m*a_x)/(sqrt(1+m*m));  // cos(atan(x)) = 1/sqrt(1+m*m), sin(pi/2-x) = cos(x)
        }

        return make_pair(alpha, rho);
    }

    T integrate_ray(
        const   diy::Master::ProxyWithLink& cp,
        const   VectorX<T>& a,
        const   VectorX<T>& b,
                bool fixed_length) const
    {
        const double pi = 3.14159265358979;
        const bool verbose = false;

        // TODO: This is for 2d only right now
        if (a.size() != 2 && b.size() != 2)
        {
            cerr << "ERROR: Incorrect dimension in integrate ray. Exiting." << endl;
            exit(1);
        }

        auto ar_coords = dualCoords(a, b);
        T alpha = ar_coords.first;
        T rho   = ar_coords.second;

        T a_x = a(0);
        T a_y = a(1);
        T b_x = b(0);
        T b_y = b(1);

        T x0, x1, y0, y1;   // end points of full line
        T u0 = 0, u1 = 0;
        T length = 0;
        T r_lim = ray_bounds_maxs(1);   // WARNING TODO: make r_lim query-able in RayMFA class
        if (fixed_length)
        {
            x0 = rho * cos(alpha) - r_lim * sin(alpha);
            x1 = rho * cos(alpha) + r_lim * sin(alpha);
            y0 = rho * sin(alpha) + r_lim * cos(alpha);
            y1 = rho * sin(alpha) - r_lim * cos(alpha);
        }
        else
        {
            get_box_intersections(alpha, rho, x0, y0, x1, y1, this->box_mins, this->box_maxs);
        }

        // parameter values along ray for 'start' and 'end'
        // compute in terms of Euclidean distance to avoid weird cases
        //   when line is nearly horizontal or vertical
        T x_sep = abs(x1 - x0);
        T y_sep = abs(y1 - y0);
        if (fixed_length)
            length = 2 * r_lim;
        else
            length = sqrt(x_sep*x_sep + y_sep*y_sep);
        
        if (x_sep > y_sep)  // want to avoid dividing by near-epsilon numbers
        {
            u0 = abs(a_x - x0) / x_sep;
            u1 = abs(b_x - x0) / x_sep;
        }
        else
        {
            u0 = abs(a_y - y0) / y_sep;
            u1 = abs(b_y - y0) / y_sep;
        }

        // Scalar valued path integrals do not have an orientation, so we always
        // want the limits of integration to go from smaller to larger.
        if (u0 > u1)
        {
            T temp  = u1;
            u1 = u0;
            u0 = temp;
        }

        if (verbose)
        {
            cerr << "RAY: (" << a(0) << ", " << a(1) << ") ---- (" << b(0) << ", " << b(1) << ")" << endl;
            cerr << "|  m: " << ((a_x==b_x) ? "inf" : to_string((b_y-a_y)/(b_x-a_x)).c_str()) << endl;
            cerr << "|  alpha:  " << alpha << ",   rho: " << rho << endl;
            cerr << "|  length: " << length << endl;
            cerr << "|  u0: " << u0 << ",  u1: " << u1 << endl;
            cerr << "+---------------------------------------\n" << endl;
        }

        VectorX<T> output(1); // todo: this is hardcoded for the first (scalar) variable only
        integrate_axis_ray(cp, alpha, rho, u0, u1, length, output);

        return output(0);
    }

    void integrate_axis_ray(
        const diy::Master::ProxyWithLink&   cp,
        T                                   alpha,
        T                                   rho,
        T                                   u0,
        T                                   u1,
        T                                   scale,
        VectorX<T>&                         output) const
    {
        T alpha_param = (alpha - ray_bounds_mins(2)) / (ray_bounds_maxs(2) - ray_bounds_mins(2));
        T rho_param = (rho - ray_bounds_mins(1)) / (ray_bounds_maxs(1) - ray_bounds_mins(1));
        
        // TODO: this is first science variable only
        ray_mfa->IntegrateAxisRay(ray_mfa->var(0), alpha_param, rho_param, u0, u1, output);

        output *= scale;
    }

    // Compute segment errors in a RayMFA
    void compute_sinogram(
                const   diy::Master::ProxyWithLink& cp,
                T extent) const
    {
        ofstream sinotruefile;
        ofstream sinoapproxfile;
        ofstream sinoerrorfile;
        string sino_true_filename = "sinogram_true_gid" + to_string(cp.gid()) + ".txt";
        string sino_approx_filename = "sinogram_approx_gid" + to_string(cp.gid()) + ".txt";
        string sino_error_filename = "sinogram_error_gid" + to_string(cp.gid()) + ".txt";
        sinotruefile.open(sino_true_filename);
        sinoapproxfile.open(sino_approx_filename);
        sinoerrorfile.open(sino_error_filename);
        int test_n_alpha = 150;
        int test_n_rho = 150;
        T r_lim = this->bounds_maxs(1);   // WARNING TODO: make r_lim query-able in RayMFA class

        int old_dom_dim = ray_dom_dim - 1;  // Assume here that dom_dim has already been incremented
        VectorX<T> start_pt(old_dom_dim), end_pt(old_dom_dim);

        for (int i = 0; i < test_n_alpha; i++)
        {
            for (int j = 0; j < test_n_rho; j++)
            {
                T alpha = 3.14159265 / (test_n_alpha-1) * i;
                T rho = r_lim*2 / (test_n_rho-1) * j - r_lim;
                T x0, x1, y0, y1;   // end points of full line

                get_box_intersections(alpha, rho, x0, y0, x1, y1, this->box_mins, this->box_maxs);
                if (x0==0 && y0==0 && x1==0 && y1==0)
                {   
                    sinotruefile << alpha << " " << rho << " " << " 0 0" << endl;
                    sinoapproxfile << alpha << " " << rho << " " << " 0 0" << endl;
                    sinoerrorfile << alpha << " " << rho << " " << " 0 0" << endl;
                }
                else
                {
                    T length = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
                    start_pt(0) = x0; 
                    start_pt(1) = y0;
                    end_pt(0) = x1;
                    end_pt(1) = y1;

                    T test_result = integrate_ray(cp, start_pt, end_pt, 1) / length;   // normalize by segment length
                    T test_actual = sintest(start_pt, end_pt) / length;

                    T e_abs = abs(test_result - test_actual);
                    T e_rel = e_abs/extent;

                    sinotruefile << alpha << " " << rho << " 0 " << test_actual << endl;
                    sinoapproxfile << alpha << " " << rho << " 0 " << test_result << endl;
                    sinoerrorfile << alpha << " " << rho << " 0 " << e_abs << endl;
                }
                
            }
        }
        sinotruefile.close();
        sinoapproxfile.close();
        sinoerrorfile.close();
        
        return;
    }

    void print_knots_ctrl(const mfa::MFA_Data<T>& model) const
    {
        VectorXi tot_nctrl_pts_dim = VectorXi::Zero(model.dom_dim);        // number contrl points per dim.
        size_t tot_nctrl_pts = 0;                                        // total number of control points

        for (auto j = 0; j < model.ntensors(); j++)
        {
            tot_nctrl_pts_dim += model.tmesh.tensor_prods[j].nctrl_pts;
            tot_nctrl_pts += model.tmesh.tensor_prods[j].nctrl_pts.prod();
        }
        // print number of control points per dimension only if there is one tensor
        if (model.ntensors() == 1)
            cerr << "# output ctrl pts     = [ " << tot_nctrl_pts_dim.transpose() << " ]" << endl;
        cerr << "tot # output ctrl pts = " << tot_nctrl_pts << endl;

        cerr << "# output knots        = [ ";
        for (auto j = 0 ; j < model.tmesh.all_knots.size(); j++)
        {
            cerr << model.tmesh.all_knots[j].size() << " ";
        }
        cerr << "]" << endl;
    }

    void print_ray_model(const diy::Master::ProxyWithLink& cp,
            bool                              error) const    // error was computed
    {
        fprintf(stderr, "gid = %d", cp.gid());
        if (!ray_mfa)
        {
            // Continue gracefully if this block never created an MFA for some reason
            fprintf(stderr, ": No Ray MFA found.\n");
            return;
        }
        else
        {
            fprintf(stderr, "\n");
        }

        // max errors over all science variables, and variables where the max error occurs
        // Initializing to 0 so we don't have to initialize with first var's error below
        T all_max_err = 0, all_max_norm_err = 0, all_max_sum_sq_err = 0, all_max_rms_err = 0, all_max_norm_rms_err = 0;
        int all_max_var = 0, all_max_norm_var = 0, all_max_sum_sq_var = 0, all_max_rms_var = 0, all_max_norm_rms_var = 0;   

        // geometry
        cerr << "\n------- geometry model -------" << endl;
        print_knots_ctrl(ray_mfa->geom());
        cerr << "-----------------------------" << endl;

        // science variables
        cerr << "\n----- science variable models -----" << endl;
        for (int i = 0; i < ray_mfa->nvars(); i++)
        {
            cerr << "\n---------- var " << i << " ----------" << endl;
            print_knots_ctrl(ray_mfa->var(i));
            cerr << "-----------------------------" << endl;

            int min_dim = ray_mfa->var(i).min_dim;
            int vardim  = ray_mfa->var_dim(i);
            MatrixX<T> varcoords = ray_input->domain.middleCols(min_dim, vardim);

            // range_extents_max is a vector containing the range extent in each component of the science variable
            // So, the size of 'range_extents_max' is the dimension of the science variable
            VectorX<T> range_extents_max = varcoords.colwise().maxCoeff() - varcoords.colwise().minCoeff();

            // 'range_extents' is the norm of the difference between the largest and smallest components
            // in each vector component. So, for each vector component we take find the difference between the
            // largest and smallest values in that component. Then we take the norm of the vector that has
            // this difference in each coordinate.
            T range_extent = range_extents_max.norm();

            if (error)
            {
                T rms_err = sqrt(ray_sum_sq_errs[i] / (ray_input->npts));
                fprintf(stderr, "range extent          = %e\n",  range_extent);
                fprintf(stderr, "max_err               = %e\n",  ray_max_errs[i]);
                fprintf(stderr, "normalized max_err    = %e\n",  ray_max_errs[i] / range_extent);
                fprintf(stderr, "sum of squared errors = %e\n",  ray_sum_sq_errs[i]);
                fprintf(stderr, "RMS error             = %e\n",  rms_err);
                fprintf(stderr, "normalized RMS error  = %e\n",  rms_err / range_extent);

                // find max over all science variables
                if (ray_max_errs[i] > all_max_err)
                {
                    all_max_err = ray_max_errs[i];
                    all_max_var = i;
                }
                if (ray_max_errs[i] / range_extent > all_max_norm_err)
                {
                    all_max_norm_err = ray_max_errs[i] / range_extent;
                    all_max_norm_var = i;
                }
                if (ray_sum_sq_errs[i] > all_max_sum_sq_err)
                {
                    all_max_sum_sq_err = ray_sum_sq_errs[i];
                    all_max_sum_sq_var = i;
                }
                if (rms_err > all_max_rms_err)
                {
                    all_max_rms_err = rms_err;
                    all_max_rms_var = i;
                }
                if (rms_err / range_extent > all_max_norm_rms_err)
                {
                    all_max_norm_rms_err = rms_err / range_extent;
                    all_max_norm_rms_var = i;
                }
            }
            cerr << "-----------------------------" << endl;
        }

        if (error)
        {
            fprintf(stderr, "\n");
            fprintf(stderr, "Maximum errors over all science variables:\n");
            fprintf(stderr, "max_err                (var %d)    = %e\n",  all_max_var,          all_max_err);
            fprintf(stderr, "normalized max_err     (var %d)    = %e\n",  all_max_norm_var,     all_max_norm_err);
            fprintf(stderr, "sum of squared errors  (var %d)    = %e\n",  all_max_sum_sq_var,   all_max_sum_sq_err);
            fprintf(stderr, "RMS error              (var %d)    = %e\n",  all_max_rms_var,      all_max_rms_err);
            fprintf(stderr, "normalized RMS error   (var %d)    = %e\n",  all_max_norm_rms_var, all_max_norm_rms_err);
        }

       cerr << "\n-----------------------------------" << endl;

        fprintf(stderr, "# input points        = %ld\n", ray_input->npts);
        fprintf(stderr, "compression ratio     = %.2f\n", this->compute_ray_compression());
    }

    // compute compression ratio
    float compute_ray_compression() const
    {
        float in_coords = (ray_input->npts) * (ray_input->pt_dim);
        float out_coords = 0.0;
        for (auto j = 0; j < ray_mfa->geom().ntensors(); j++)
            out_coords += ray_mfa->geom().tmesh.tensor_prods[j].ctrl_pts.rows() *
                ray_mfa->geom().tmesh.tensor_prods[j].ctrl_pts.cols();
        for (auto j = 0; j < ray_mfa->geom().tmesh.all_knots.size(); j++)
            out_coords += ray_mfa->geom().tmesh.all_knots[j].size();
        for (auto i = 0; i < ray_mfa->nvars(); i++)
        {
            for (auto j = 0; j < ray_mfa->var(i).ntensors(); j++)
                out_coords += ray_mfa->var(i).tmesh.tensor_prods[j].ctrl_pts.rows() *
                    ray_mfa->var(i).tmesh.tensor_prods[j].ctrl_pts.cols();
            for (auto j = 0; j < ray_mfa->var(i).tmesh.all_knots.size(); j++)
                out_coords += ray_mfa->var(i).tmesh.all_knots[j].size();
        }
        return in_coords / out_coords;
    }
}; // RayBlock

#endif // _MFA_RAY_BLOCK_HPP