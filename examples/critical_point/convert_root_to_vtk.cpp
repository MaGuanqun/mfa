//-----------
//Compute derivative control points of an MFA for an entrie block
#include "mfa/mfa.hpp"

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

#include "block.hpp"

// #include "EigenMatrixIO.h"
#include "convert/writer.hpp"
// #include "save_control_data.hpp"
#include "utility_function.h"
#include"find_unique_root.h"


template<typename T>
bool testInBlock(VectorX<T>& point, VectorX<T>& min, VectorX<T>& max)
{
    for(int i=0;i<min.size();++i)
    {
        if(point(i)<min(i) || point(i)>max(i))
        {
            return false;
        }
    }
    return true;
}


template<typename T>
void proceed_find_unique_root(std::vector<Eigen::VectorX<T>>& points, VectorX<T>& min, VectorX<T>& max,std::vector<Eigen::VectorX<T>>& unique_points, T same_root_epsilon)
{

       // normalize point
    // for (auto& point : points) {
    //     point(0) = (point(0) - min(0)) / (max(0) - min(0));
    //     point(1) = (point(1) - min(1)) / (max(1) - min(1));
    // }

    // Find the unique points

    
    spatial_hashing::find_all_unique_root(points,unique_points,same_root_epsilon);

    // Write the unique points to a file
    std::cout<<"unique_points.size() "<<unique_points.size()<<std::endl;
    std::cout<<"points.size() "<<points.size()<<std::endl;
    std::cout<<"same root epsilon is: "<<same_root_epsilon<<std::endl;


    // for (auto& point : unique_points) {
    //     point(0) = point(0) * (max(0) - min(0))+min(0);
    //     point(1) = point(1) * (max(1) - min(1))+min(1);
    // }
}



int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "derivative_root_point.dat";               // diy input file
    string outfile = "root.csv";

    string duplicate_file;

    // default command line arguments
    int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help


    string function_value_input;

    string input_mfa_file = "approx.mfa";
    
    string input_shrink_ratio = "0-1-0-1-0-1";

    double same_root_epsilon = SAME_ROOT_EPSILON;

    int reduce_duplication = 0;
    int input_vector_vector = 0;

    // get command line arguments
    opts::Options ops;
    //ops >> opts::Option('d', "deriv",   deriv,   " which derivative to take (1 = 1st, 2 = 2nd, ...)");
    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    //ops >> opts::Option('a', "partial", partial, " dimension of 1 partial derivative only");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('u', "duplicate_file",    duplicate_file,    " deplicate degree of root");
    ops >> opts::Option('v', "function_value_input",    function_value_input,    " set a constant function value");
    ops >> opts::Option('o', "outfile", outfile,    " current the csv file");

    ops >> opts::Option('s', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-z1-z2-...\"");
    ops >> opts::Option('i', "input_mfa_file",      input_mfa_file,     " input mfa file name");
    ops >> opts::Option('e', "same_root_epsilon",      same_root_epsilon,     " same root epsilon");
    ops >> opts::Option('d', "reduce_duplication",      reduce_duplication,     " 0 is not reducing duplication, 1 is reducing duplication");
    ops >> opts::Option('p', "vector_vector",    input_vector_vector,       " if the input data need to be transferred to vector<vector<>>");

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
    
    std::vector<Eigen::MatrixXd> root;

    std::ifstream file(infile.c_str());
    if (!file) {
        std::cerr << "File does not exist: "<< std::endl;
        return 0;
    }


    utility::loadMatrixVector(infile.c_str(),root);

    std::vector<Eigen::MatrixXd> duplicate_num;
    if(!duplicate_file.empty())
    {
        utility::loadMatrixVector(duplicate_file.c_str(),duplicate_num);
    }

    std::cout<<"root size is: "<<root.size()<<std::endl;
   

    std::istringstream iss(input_shrink_ratio);

    std::cout<<"shrink range is: "<<input_shrink_ratio<<std::endl;

    std::vector<double> shrink_ratio;
    double number;
    std::string token;
    while (std::getline(iss, token, '-')) {
        std::istringstream tokenStream(token);
        if (tokenStream >> number) {
            shrink_ratio.push_back(number);
        }
    } 

    std::cout<<"shrink ratio is: ";
    for(auto num:shrink_ratio)
    {
        std::cout<<num<<" ";
    }
    std::cout<<std::endl;

    std::cout<<"threshold "<<same_root_epsilon<<std::endl;

   // initialize DIY
    diy::FileStorage storage("./DIY.XXXXXX");     // used for blocks to be moved out of core
    diy::Master      master(world,
            1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy);
    diy::ContiguousAssigner   assigner(world.size(), -1); // number of blocks set by read_blocks()

    diy::io::read_blocks(input_mfa_file.c_str(), world, assigner, master, &Block<real_t>::load);
    std::cout << master.size() << " blocks read from file "<< input_mfa_file << "\n\n";

    

    std::vector<std::vector<Eigen::VectorXd>> filterd_root(root.size());    

    VectorXd min,max;
    VectorXd min_ori,max_ori;
    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {
        min.resize(b->core_mins.size());
        max.resize(b->core_mins.size());
        for(int i=0;i<b->core_mins.size();++i)
        {
            if(i>=shrink_ratio.size()/2)
            {
                shrink_ratio.push_back(0);
                shrink_ratio.push_back(1);
            }

            shrink_ratio[2*i] = round((b->input->ndom_pts(i)-1)*shrink_ratio[2*i])/(b->input->ndom_pts(i)-1);
            shrink_ratio[2*i+1] = round((b->input->ndom_pts(i)-1)*shrink_ratio[2*i+1])/(b->input->ndom_pts(i)-1);

            min(i) = b->core_mins(i) + shrink_ratio[2*i]*(b->core_maxs(i)-b->core_mins(i));
            max(i) = b->core_mins(i) + shrink_ratio[2*i+1]*(b->core_maxs(i)-b->core_mins(i));

            min_ori=b->core_mins;
            max_ori=b->core_maxs;

        }
         std::cout<<" mfa min max are: "<<b->core_mins.transpose()<<" "<<b->core_maxs.transpose()<<std::endl;
    });

   

    std::cout<<"restored min max are: "<<min.transpose()<<" "<<max.transpose()<<std::endl;

    int start_place=0;
    if(input_vector_vector)
    {
        start_place=root[0].data()[0]+1;
    }
    std::cout<<"initial root size "<<root[0].rows()<<" "<<root[0].cols()<<std::endl;

    for(int i=0;i<root.size();++i)
    {
        filterd_root[i].reserve(root[i].rows());
        for(int j=start_place;j<root[i].rows();++j)
        {
            Eigen::VectorXd temp;
            if(function_value_input.empty())
                temp = root[i].row(j).transpose();//.cast<float>().transpose();
            else
                temp = root[i].row(j).head(root[i].cols()-1).transpose();//.cast<float>().transpose();
            // if(testInBlock(temp,min,max))
            // {
                filterd_root[i].push_back(temp);
            // }
        }
    }


    std::vector<std::vector<Eigen::VectorXd>> unique_root(filterd_root.size());

    if(reduce_duplication==1)
    {
        proceed_find_unique_root(filterd_root[0],min,max,unique_root[0],same_root_epsilon);
    }
    else{
        unique_root = filterd_root;
    }
    // 
    

    std::cout<<"filterd root size "<<unique_root[0].size()<<std::endl;


    if(!function_value_input.empty())
    {
        utility::saveToCSV(outfile, unique_root[0],std::stod(function_value_input));
    }
    else
    {
        utility::saveToCSV(outfile, unique_root[0]);
    }
    


    // std::vector<Eigen::MatrixXf> combine_root(root.size());
    // for(int i=0;i<root.size();++i)
    // {

    //     if(root[i].cols()<4)
    //     {
    //         combine_root[i].resize(root[i].cols(),unique_root[i].size());            
    //         //set the position to be value
    //         for(int j=0;j<unique_root[i].size();++j)
    //         {
    //             combine_root[i].col(j) = unique_root[i][j].cast<float>();
    //         }

    //         // memcpy(combine_root[i].data(),unique_root[i][0].data(),sizeof(float)*combine_root[i].rows()*combine_root[i].cols());

    //         if (!function_value_input.empty())
    //         {
    //             combine_root[i].row(combine_root[i].rows()-1).setConstant(std::stof(function_value_input));
    //         }

    //         //set the position to be zero
    //         // combine_root[i].block(0,0,root[i].cols()-1,root[i].rows()) = root[i].block(0,0,root[i].rows(),root[i].cols()-1).transpose().cast<float>();  
    //     }
    //     else
    //     {
    //         combine_root[i].resize(3,unique_root[i].size());            
    //         combine_root[i].setZero();

    //         for(int j=0;j<unique_root[i].size();++j)
    //         {
    //             combine_root[i].col(j) = unique_root[i][j].head(3).cast<float>();
    //         }

    //         // combine_root[i].block(0,0,3,root[i].rows())= root[i].block(0,0,root[i].rows(),3).transpose().cast<float>();

    //        // combine_root[i].block(2,0,1,root[i].rows())= root[i].col(root[i].cols()-1).transpose().cast<float>();
    //     }  
    // }





    // int nvars=1;
    // std::vector<int>vardims(nvars);
    // std::vector<int>centerings(nvars);
    // char** varnames     = new char*[nvars];
    // float** pt_data = new float*[nvars];

    // for (int i = 0; i < nvars; i++)
    // {
    //     vardims[i]=1;
    //     centerings[i]   = 1;
    //     pt_data[i] = new float[combine_root[i].cols()];
    //     memset(pt_data[i],0,sizeof(float)*combine_root[i].cols()); //this should be the value of every variable, they are root here set them to be zero

    //     varnames[i]     = new char[256];
    //     sprintf(varnames[i], "var%d", i);

    //     if(duplicate_file.empty())
    //     {
    //         for(int j=0;j<combine_root[i].cols();j++)
    //         {
    //             pt_data[i][j] = combine_root[i](combine_root[i].rows()-1,j);
    //         }            
    //     }
    //     else
    //     {
    //         for(int j=0;j<combine_root[i].cols();j++)
    //         {
    //             pt_data[i][j] = duplicate_num[i](j,0);
    //         }     
    //     }

    // }




    // // for(int i=0;i<root[0].cols();i++)
    // // {
    // //     memcpy(pt_data,root[0].data()+root[0].cols()*root[0].rows(),sizeof(float) * control_point[i][m].rows());
    // // }


    // // std::cout<<combine_root[0]<<std::endl;

    // write_point_mesh(outfile.c_str(),0,combine_root[0].cols(),combine_root[0].data(),
    // nvars,vardims.data(),varnames,pt_data);//


}




