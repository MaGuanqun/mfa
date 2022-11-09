//--------------------------------------------------------------
// example of encoding / decoding higher dimensional data w/ fixed number of control points and a
// single block in a split model w/ one model containing geometry and other model science variables
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <set>

#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include "opts.h"

#include "block.hpp"
#include "example-setup.hpp"

using namespace std;

int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    int nblocks     = 1;                        // number of local blocks
    int tot_blocks  = nblocks * world.size();   // number of global blocks
    int mem_blocks  = -1;                       // everything in core for now
    int num_threads = 1;                        // needed in order to do timing

    // default command line arguments
    int         pt_dim          = 3;        // dimension of input points
    int         dom_dim         = 2;        // dimension of domain (<= pt_dim)
    int         scalar          = 1;        // flag for scalar or vector-valued science variables (0 == multiple scalar vars)
    int         geom_degree     = 1;        // degree for geometry (same for all dims)
    int         vars_degree     = 4;        // degree for science variables (same for all dims)
    int         ndomp           = 100;      // input number of domain points (same for all dims)
    int         ntest           = 0;        // number of input test points in each dim for analytical error tests
    int         geom_nctrl      = -1;       // input number of control points for geometry (same for all dims)
    vector<int> vars_nctrl      = {11};     // initial # control points for all science variables (default same for all dims)
    string      input           = "sinc";   // input dataset
    int         weighted        = 1;        // solve for and use weights (bool 0/1)
    real_t      rot             = 0.0;      // rotation angle in degrees
    real_t      twist           = 0.0;      // twist (waviness) of domain (0.0-1.0)
    real_t      noise           = 0.0;      // fraction of noise
    int         error           = 1;        // decode all input points and check error (bool 0/1)
    string      infile;                     // input file name
    string      infile2;
    int         structured      = 1;        // input data format (bool 0/1)
    int         rand_seed       = -1;       // seed to use for random data generation (-1 == no randomization)
    real_t      regularization  = 0;        // smoothing parameter for models with non-uniform input density (0 == no smoothing)
    int         reg1and2        = 0;        // flag for regularizer: 0 = regularize only 2nd derivs. 1 = regularize 1st and 2nd
    int         verbose         = 1;        // MFA verbosity (0 = no extra output)
    vector<int> decode_grid     = {};       // Grid size for uniform decoding
    bool        help            = false;    // show help

    // Constants for this example
    const bool adaptive = false;

    // get command line arguments
    opts::Options ops;
    ops >> opts::Option('d', "pt_dim",      pt_dim,     " dimension of points");
    ops >> opts::Option('m', "dom_dim",     dom_dim,    " dimension of domain");
    ops >> opts::Option('l', "scalar",      scalar,     " flag for scalar or vector-valued science variables");
    ops >> opts::Option('p', "geom_degree", geom_degree," degree in each dimension of geometry");
    ops >> opts::Option('q', "vars_degree", vars_degree," degree in each dimension of science variables");
    ops >> opts::Option('n', "ndomp",       ndomp,      " number of input points in each dimension of domain");
    ops >> opts::Option('a', "ntest",       ntest,      " number of test points in each dimension of domain (for analytical error calculation)");
    ops >> opts::Option('g', "geom_nctrl",  geom_nctrl, " number of control points in each dimension of geometry");
    ops >> opts::Option('v', "vars_nctrl",  vars_nctrl, " number of control points in each dimension of all science variables");
    ops >> opts::Option('i', "input",       input,      " input dataset");
    ops >> opts::Option('w', "weights",     weighted,   " solve for and use weights");
    ops >> opts::Option('r', "rotate",      rot,        " rotation angle of domain in degrees");
    ops >> opts::Option('t', "twist",       twist,      " twist (waviness) of domain (0.0-1.0)");
    ops >> opts::Option('s', "noise",       noise,      " fraction of noise (0.0 - 1.0)");
    ops >> opts::Option('c', "error",       error,      " decode entire error field (default=true)");
    ops >> opts::Option('f', "infile",      infile,     " input file name");
    ops >> opts::Option('h', "help",        help,       " show help");
    ops >> opts::Option('x', "structured",  structured, " input data format (default=structured=true)");
    ops >> opts::Option('y', "rand_seed",   rand_seed,  " seed for random point generation (-1 = no randomization, default)");
    ops >> opts::Option('b', "regularization", regularization, "smoothing parameter for models with non-uniform input density");
    ops >> opts::Option('k', "reg1and2",    reg1and2,   " regularize both 1st and 2nd derivatives (if =1) or just 2nd (if =0)");
    ops >> opts::Option('u', "grid_decode", decode_grid," size of regular grid to decode MFA");
    ops >> opts::Option('z', "infile2",     infile2,    " extra data file (some apps require two file paths");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // print input arguments
    echo_mfa_settings("fixed example", pt_dim, dom_dim, scalar, geom_degree, geom_nctrl, vars_degree, vars_nctrl,
                        regularization, reg1and2, weighted, adaptive, 0, 0);
    echo_data_settings(ndomp, ntest, input, infile, noise, rot, twist, structured, rand_seed);

    // initialize DIY
    diy::FileStorage          storage("./DIY.XXXXXX"); // used for blocks to be moved out of core
    diy::Master               master(world,
                                     num_threads,
                                     mem_blocks,
                                     &Block<real_t>::create,
                                     &Block<real_t>::destroy,
                                     &storage,
                                     &Block<real_t>::save,
                                     &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), tot_blocks);

    // set global domain bounds and decompose
    Bounds<real_t> dom_bounds(dom_dim);
    set_dom_bounds(dom_bounds, input);
    
    Decomposer<real_t> decomposer(dom_dim, dom_bounds, tot_blocks);
    decomposer.decompose(world.rank(),
                         assigner,
                         [&](int gid, const Bounds<real_t>& core, const Bounds<real_t>& bounds, const Bounds<real_t>& domain, const RCLink<real_t>& link)
                         { Block<real_t>::add(gid, core, bounds, domain, link, master, dom_dim, pt_dim, 0.0); });

    // If scalar == true, assume all science vars are scalar. Else one vector-valued var
    // We assume that dom_dim == geom_dim
    // Different examples can reset this below
    vector<int> model_dims;
    if (scalar) // Set up (pt_dim - dom_dim) separate scalar variables
    {
        model_dims.assign(pt_dim - dom_dim + 1, 1);
        model_dims[0] = dom_dim;                        // index 0 == geometry
    }
    else    // Set up a single vector-valued variable
    {   
        model_dims = {dom_dim, pt_dim - dom_dim};
    }

    // Create empty info classes
    MFAInfo     mfa_info(dom_dim, verbose);
    DomainArgs  d_args(dom_dim, model_dims);
    
    // set up parameters for examples
    setup_args(dom_dim, pt_dim, model_dims, geom_degree, geom_nctrl, vars_degree, vars_nctrl,
                input, infile, infile2, ndomp, structured, rand_seed, rot, twist, noise,
                weighted, reg1and2, regularization, adaptive, verbose, mfa_info, d_args);

    // Create data set for modeling. Input keywords are defined in example-setup.hpp
    if (analytical_signals.count(input) == 1)
    {
        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        { 
            b->generate_analytical_data(cp, input, mfa_info, d_args); 
        });
    }
    else if (datasets_3d.count(input) == 1)
    {
        if (dom_dim > 3)
        {
            fprintf(stderr, "\'%s\' data only available with dimension <= 3\n", input);
            exit(0);
        }

        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        { 
            if (dom_dim == 1) b->read_1d_slice_3d_vector_data(cp, mfa_info, d_args);
            if (dom_dim == 2) b->read_2d_slice_3d_vector_data(cp, mfa_info, d_args);
            if (dom_dim == 3) b->read_3d_vector_data(cp, mfa_info, d_args);
        });
        // for testing, hard-code a subset of a 3d domain, 1/2 the size in each dim and centered
        // master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        //         { b->read_3d_subset_3d_vector_data(cp, d_args); });
    }
    else if (datasets_2d.count(input) == 1)
    {
        if (dom_dim != 2)
        {
            fprintf(stderr, "\'%s\' data only available with dimension 2\n");
            exit(0);
        }

        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        {
            b->read_2d_scalar_data(cp, mfa_info, d_args); 
        });
    }
    else if (datasets_unstructured.count(input) == 1)
    {
        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        {
            b->read_3d_unstructured_data(cp, mfa_info, d_args); 
        });
    }
    else if (input=="xgc")
    {
        string xgc_coords_name = "coordinates/values";
        string xgc_values_name = "dneOverne0";
        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        {
            b->read_2d_hdf5_data<float>(cp, mfa_info, d_args, xgc_coords_name, xgc_values_name);
        });
    }
    else
    {
        cerr << "Input keyword \'" << input << "\' not recognized. Exiting." << endl;
        exit(0);
    }

    // compute the MFA
    fprintf(stderr, "\nStarting fixed encoding...\n\n");
    double encode_time = MPI_Wtime();
    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
            { b->fixed_encode_block(cp, mfa_info); });
    encode_time = MPI_Wtime() - encode_time;
    fprintf(stderr, "\n\nFixed encoding done.\n\n");

    // debug: compute error field for visualization and max error to verify that it is below the threshold
    double decode_time = MPI_Wtime();
    if (error)
    {
        fprintf(stderr, "\nFinal decoding and computing max. error...\n");
#ifdef CURVE_PARAMS     // normal distance
        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                { b->error(cp, 1, true); });
#else                   // range coordinate difference
        bool saved_basis = structured; // TODO: basis functions are currently only saved during encoding of structured data
        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        { 
            b->range_error(cp, true, saved_basis);
        });
#endif
        decode_time = MPI_Wtime() - decode_time;
    }
    else if (decode_grid.size() == dom_dim)
    {
        fprintf(stderr, "\nDecoding on regular grid of size %s\n", mfa::print_vec(decode_grid).c_str());
        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        {
            b->decode_block_grid(cp, decode_grid);
        });
        decode_time = MPI_Wtime() - decode_time;
    }

    // debug: write original and approximated data for reading into z-checker
    // only for one block (one file name used, ie, last block will overwrite earlier ones)
//     master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
//             { b->write_raw(cp); });

    // debug: save knot span domains for comparing error with location in knot span
//     master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
//             { b->knot_span_domains(cp); });

    // compute the norms of analytical errors synthetic function w/o noise at different domain points than the input
    if (ntest > 0)
    {
        cerr << "Computing analytical error" << endl;
        int nvars = model_dims.size() - 1;
        vector<real_t> L1(nvars), L2(nvars), Linf(nvars);                                // L-1, 2, infinity norms
        d_args.ndom_pts = vector<int>(dom_dim, ntest);
        mfa::PointSet<real_t>* temp_in = nullptr;
        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                { b->analytical_error_field(cp, input, L1, L2, Linf, d_args,  temp_in, b->approx, b->errs); });

        // print analytical errors
        for (int i = 0; i < nvars; i++)
        {
            fprintf(stderr, "\n------ Analytical error: Var %i -------\n", i);
            fprintf(stderr, "L-1        norm = %e\n", L1[i]);
            fprintf(stderr, "L-2        norm = %e\n", L2[i]);
            fprintf(stderr, "L-infinity norm = %e\n", Linf[i]);
            fprintf(stderr, "-------------------------------------\n\n");
        }

    }

    // print results
    fprintf(stderr, "\n------- Final block results --------\n");
    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
            { b->print_block(cp, error); });
    fprintf(stderr, "encoding time         = %.3lf s.\n", encode_time);
    if (error)
        fprintf(stderr, "decoding time         = %.3lf s.\n", decode_time);
    fprintf(stderr, "-------------------------------------\n\n");

    // save the results in diy format
    diy::io::write_blocks("approx.mfa", world, master);
}
