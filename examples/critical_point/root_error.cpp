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

// #include "EigenMatrixIO.h"
#include "convert/writer.hpp"
// #include "save_control_data.hpp"
#include "utility_function.h"
#include"find_all_root.h"

using namespace std;

int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file
    string inroot = "root.dat";
    // string outfile = "root.vtk";

    // default command line arguments
    int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help

    // get command line arguments
    opts::Options ops;
    //ops >> opts::Option('d', "deriv",   deriv,   " which derivative to take (1 = 1st, 2 = 2nd, ...)");
    ops >> opts::Option('f', "infile",  infile,  " diy input approx MFA file");
    ops >> opts::Option('i', "root",  inroot,  " diy input root file");
    //ops >> opts::Option('a', "partial", partial, " dimension of 1 partial derivative only");
    ops >> opts::Option('h', "help",    help,    " show help");
   // ops >> opts::Option('o', "outfile", outfile,    " control points of first derivative");

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
    
    std::vector<Eigen::MatrixXd> ori_root;




    utility::loadMatrixVector(inroot.c_str(),ori_root);

    std::vector<std::vector<Eigen::VectorXd>> root;
    root.resize(ori_root.size());
    for(int i=0;i<root.size();++i)
    {
        int domain_num = ori_root[i].cols()-1;
        root[i].resize(ori_root[i].rows());
        for(int j=0;j<root[i].size();j++)
        {
            root[i][j]=ori_root[i].block(j,0,1,domain_num).transpose();
        }
    }

  
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
    int index = 0;

    std::vector<std::vector<VectorX<double>>> value(master.size());//[blocks,...]

    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {

        find_all_roots::getDerivative(b->mfa,b,root[index],value[index]);

        index++;    
    });

    for(int i=0;i<master.size();i++)
    {
        for(int j=0;j<value[i].size();j++)
        {
            if((root[i][j][0]<0.01 || root[i][j][0]>0.99)||(root[i][j][1]<0.01 || root[i][j][1]>0.99))
            {
                std::cout<<root[i][j].transpose()<<" "<<value[i][j].transpose()<<std::endl;
            }            
        }
    }


}




