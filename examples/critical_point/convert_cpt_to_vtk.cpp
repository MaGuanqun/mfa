//-----------
//Compute derivative control points of an MFA for an entrie block
#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

#include <diy/master.hpp>
// #include <diy/reduce-operations.hpp>
// #include <diy/decomposition.hpp>
// #include <diy/assigner.hpp>
#include <diy/io/block.hpp>


#include "opts.h"

// #include "block.hpp"

// #include "EigenMatrixIO.h"
#include "convert/writer.hpp"
// #include ".writer.hpp"
#include "save_control_data.hpp"
#include "utility_function.h"




int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "derivative_control_point.dat";               // diy input file
    string outfile = "partial_derivative_control_point";

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
    std::vector<std::vector<std::vector<std::vector<double>>>> geo_control_point; //[deriv][vars][dom][...]
    std::vector<std::vector<MatrixX<double>>> sci_deriv_control_points;//[vars][partial_deriv][...]

    save_control_points::load_control_points(infile.c_str(),geo_control_point,sci_deriv_control_points);
    
    std::vector<std::vector<Eigen::MatrixXf>> control_point; //[vars][dom]
    control_point.resize(sci_deriv_control_points.size());

    Eigen::VectorXi ori_cpt_size,current_cpt_size;
    VectorXi number_in_every_domain;
    int domain_size;

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

    // for(int i=0;i<sci_deriv_control_points.size();++i)
    // {
    //     for(int j=0;j<sci_deriv_control_points[i].size();j++)
    //     {
    //       std::cout<<sci_deriv_control_points[i][j].transpose()<<std::endl;
    //     }
    // }


    VectorXi domain_index;
    for(int i=0;i<sci_deriv_control_points.size();i++) //vars
    {
        domain_size = geo_control_point[0][i].size();
        domain_index.resize(domain_size);
        ori_cpt_size.resize(domain_size);
        for(int j=0;j<ori_cpt_size.size();j++)
        {
            ori_cpt_size[j]=geo_control_point[0][i][j].size();
        }
        control_point[i].resize(domain_size);
        for(int j=0;j<domain_size;j++) //dom, we need to compute f_x, f_y, f_z...
        {
            current_cpt_size = ori_cpt_size;
            current_cpt_size[j]-=1;

            utility::obtain_number_in_every_domain(current_cpt_size,number_in_every_domain);

            // number_in_every_domain.resize(domain_size);
            // number_in_every_domain[0]=1;
            // for(int k=1;k<domain_size;k++)
            // {
            //     number_in_every_domain[k]=number_in_every_domain[k-1]*current_cpt_size[k-1];
            // }

            size_t pt_number = current_cpt_size.prod();

            if(domain_size<2)
            {
                control_point[i][j].resize(pt_number,3);
                control_point[i][j].setZero();
            }
            else
            {
                control_point[i][j].resize(pt_number,domain_size+1);
            }

            for(auto k=0;k<pt_number;k++)
            {
                utility::obtainDomainIndex(k,domain_index,number_in_every_domain);
                for(int m=0;m<domain_size;m++)
                {
                    //convert k to index in every dimension [x,y,...]
                    if(m==j){
                        control_point[i][j](k,m)=geo_control_point[1][i][m][domain_index[m]];
                    }
                    else
                    {
                        control_point[i][j](k,m)=geo_control_point[0][i][m][domain_index[m]];
                    }
                }
                
            }

            control_point[i][j].col(domain_size) = sci_deriv_control_points[i][j].col(0).cast<float>();
        }    
    }
    // for(int i=0;i<control_point.size();++i)
    // {
    //     for(int j=0;j<control_point[i].size();j++)
    //     {
    //       std::cout<<control_point[i][j]<<std::endl;
    //       std::cout<<"===="<<std::endl;
    //     }
    // }


    std::cout<<"derivative-control-point "<<std::endl;
    std::vector<std::vector<Eigen::MatrixXf>> transpose_control_point(control_point.size());
    for(int i=0;i<control_point.size();i++)
    {
        transpose_control_point[i].resize(control_point[i].size());
        for(int j=0;j<control_point[i].size();j++)
        {
            transpose_control_point[i][j] = control_point[i][j].transpose();

        }

    }


    domain_size = geo_control_point[0][0].size();

    std::vector<std::string>out_file_name(domain_size);

    for(int m=0;m< domain_size;m++) //domain
    {    
        int nvars=1;
        std::vector<int>vardims(nvars);
        std::vector<int>centerings(nvars);
        char** varnames     = new char*[nvars];
        float** pt_data = new float*[nvars];

        out_file_name[m] = outfile +"_" + std::to_string(m)+".vtk";

        for (int i = 0; i < nvars; i++)
        {
            vardims[i]      = 1;  
            varnames[i]     = new char[256];
            centerings[i]   = 1;
            sprintf(varnames[i], "var%d", i);
            pt_data[i]  = new float[control_point[i][m].rows()];
            memcpy(pt_data[i],control_point[i][m].data()+domain_size*control_point[i][m].rows(),sizeof(float) * control_point[i][m].rows());
        }    

        write_point_mesh(out_file_name[m].c_str(),0,control_point[0][m].rows(),transpose_control_point[0][m].data(),
        nvars,vardims.data(),varnames,pt_data);//

        for (int i = 0; i < nvars; i++)
            delete[] varnames[i];
        delete[] varnames;
        for (int j = 0; j < nvars; j++)
        {
            delete[] pt_data[j];
        }
        delete[] pt_data;
    }

}




