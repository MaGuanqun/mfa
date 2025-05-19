//-----------
//Compute derivative control points of an MFA for an entrie block
#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>



#include "opts.h"

#include "block.hpp"

#include "diff_control_point.hpp"
#include "save_control_data.hpp"
using namespace std;


int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file
    string outfile = "derivative_control_point.dat";

    // default command line arguments
    int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help

    // get command line arguments
    opts::Options ops;
    ops >> opts::Option('d', "deriv",   deriv,   " which derivative to take (1 = 1st, 2 = 2nd, ...)");
    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('a', "partial", partial, " dimension of 1 partial derivative only");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('o', "outfile", outfile,    " control points of first derivative");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // echo args
    fprintf(stderr, "\n--------- Input arguments ----------\n");
    cerr <<
        "deriv = "    << deriv << endl;
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

    
    // initialize DIY
    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks to be moved out of core
    diy::Master      master(world,
            -1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy,
            &storage,
            &Block<real_t>::save,
            &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), -1);   // number of blocks set by read_blocks()

    // read MFA model
    diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
    int nblocks = master.size();
    std::cout << nblocks << " blocks read from file "<< infile << "\n";

    //compute control point
    double decode_time = MPI_Wtime();

    std::vector<std::vector<std::vector<std::vector<real_t>>>> geo_control_point;
    std::vector<std::vector<MatrixX<real_t>>> sci_deriv_control_points;

    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
            {
                diff_control_point_block(b,deriv,geo_control_point,sci_deriv_control_points); });
    decode_time = MPI_Wtime() - decode_time;

    std::cout<<"variable number "<<geo_control_point.size()<<std::endl;

    // for(int i=0;i<geo_control_point.size();++i)
    // {
    //     for(int j=0;j<geo_control_point[i][0].size();j++)
    //     {
    //         for(int k=0;k<geo_control_point[i][0][j].size();++k)
    //         {
    //             std::cout<<geo_control_point[i][0][j][k]<<" ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }

    save_control_points::save_control_points(outfile.c_str(),geo_control_point,sci_deriv_control_points);

    // EigenMatrixIO::write_binary(outfile.c_str(),deriv_control_points[0]);



}

