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

#include <chrono>
#include    "opts.h"


#include <tbb/tbb.h>
#include    "block.hpp"


// #include "transfer_data.h"
// #include <fstream>
#include <iostream>
// #include <Eigen/Dense>
#include "ridge_valley_graph.h"
#include "find_jacobi_set.h"
#include "mfa_extend.h"

using namespace std;


template<typename T>
void fixPoint(std::vector<VectorX<T>>& point, const VectorX<T>& domain_min, const VectorX<T>& domain_max)
{
    for(auto i=0;i<point.size();++i)
    {
        for(auto j=0;j<point[i].size();++j)
        {
            if(point[i][j]<domain_min[j])
            {
                point[i][j]=domain_min[j];
            }
            if(point[i][j]>domain_max[j])
            {
                point[i][j]=domain_max[j];
            }
        }
    }
}

template<typename T>
void getError(const Block<T>* b, std::vector<VectorX<T>>& point,std::vector<T>& error,
    T theoritical_value)
{
    tbb::affinity_partitioner ap;
    fixPoint(point,b->core_mins,b->core_maxs);
    error.resize(point.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, point.size()),
    [&](const tbb::blocked_range<size_t>& range)
    {
        VectorX<T> value0(1); 
        for(auto i=range.begin();i!=range.end();++i)
        {
            mfa_extend::recover_mfa(b,point[i],value0);
            // VectorX<T> param = (point[i]-domain_min).cwiseQuotient(domain_range);
            // mfa->DecodePt(*mfa_data,param,value0);
            error[i]=std::abs(value0[0]-theoritical_value);
        }
     },ap
    );     
    // for(auto i=0;i<point.size();++i)
    // {
    //     VectorX<T> param = (point[i]-domain_min).cwiseQuotient(domain_range);
    //     mfa->DecodePt(*mfa_data,param,value0);
    //     error[i]=std::abs(value0[0]-theoritical_value);
    // }
    // for(auto i=0;i<error.size();++i)
    // {
    //     std::cout<<error[i]<<" ";
    // }

}


template<typename T>
void get_rv_error(const Block<T>* b, std::vector<VectorX<T>>& point,std::vector<T>& error)
{
    tbb::affinity_partitioner ap;
    fixPoint(point,b->core_mins,b->core_maxs);
    error.resize(point.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, point.size()),
    [&](const tbb::blocked_range<size_t>& range)
    {

        for(auto i=range.begin();i!=range.end();++i)
        {
            ridge_valley_graph::compute_h(point[i],b,error[i]);
            error[i]=std::abs(error[i]);
        }
     },ap
    );     
    // for(auto i=0;i<point.size();++i)
    // {
    //     VectorX<T> param = (point[i]-domain_min).cwiseQuotient(domain_range);
    //     mfa->DecodePt(*mfa_data,param,value0);
    //     error[i]=std::abs(value0[0]-theoritical_value);
    // }
    // for(auto i=0;i<error.size();++i)
    // {
    //     std::cout<<error[i]<<" ";
    // }

}

template<typename T>
void get_js_error(const Block<T>* b,const Block<T>* b2,  std::vector<VectorX<T>>& point,std::vector<T>& error)
{
    tbb::affinity_partitioner ap;


    fixPoint(point,b->core_mins,b->core_maxs);

    error.resize(point.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, point.size()),
    [&](const tbb::blocked_range<size_t>& range)
    {

        for(auto i=range.begin();i!=range.end();++i)
        {
            jacobi_set::compute_h(point[i],b,b2,error[i]);
            error[i]=std::abs(error[i]);
        }
     },ap
    );     



}




template<typename T>
void read_bin(string& point_file_name,std::vector<VectorX<T>>& data)
{

    // Open the binary file for reading
    std::ifstream file(point_file_name.c_str(), std::ios::binary);




    if (file) {
        // Get the size of the file
        file.seekg(0, std::ios::end);
        std::streampos fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        // Read the binary data into a buffer
        std::vector<char> buffer(fileSize);
        file.read(buffer.data(), fileSize);

        data.resize(fileSize / (3 * sizeof(T)));



        printf("data size %d\n",int(data.size()));

        for (size_t i = 0; i < data.size(); i++) {
            data[i].resize(2);
            data[i] = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<T*>(buffer.data() + i * 3 * sizeof(T)), 2);
        }
        // Clean up
        file.close();

        std::cout << "Read " << data.size() << " rows from the binary file." << std::endl;

        // for(int i=0;i<data.size();++i)
        // {
        //     std::cout<<data[i].transpose()<<std::endl;
        // }
    } else {
        std::cout << "Failed to open the binary file." << std::endl;
    }

}


int main(int argc, char** argv)
{

    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file

    string root_file = "root_point.bin";

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;

   

    string infile2 = "approx2.mfa";               // diy input file

    
    double point_func_value = 0.0; //value of the point need to be extracted.
    // double same_root_epsilon = SAME_ROOT_EPSILON;

    bool extract_isocontour = false;
    bool extract_ridge_valley = true;

    double root_finding_epsilon = 1e-8;

    double threshold_distance_stop_tracing = 0.3;
    string input_shrink_ratio = "0-1-0-1";

    string func_error_file = "func_error.dat";

    int func_type = 0;

    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('g', "infile2",  infile2,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('r', "root_file", root_file,    " root");
    ops >> opts::Option('s', "func_error_file",    func_error_file,       " file name of function error");
    ops >> opts::Option('v', "point_func_value",    point_func_value,       " value of the point need to be extracted");
    ops >> opts::Option('z', "func_type",    func_type,       " 0: contour, 1: ridge valley graph, 2: Jacobi set");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }





    // initialize DIY
    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks to be moved out of core
    diy::Master      master(world,
            -1,
            -1,
            &Block<double>::create,
            &Block<double>::destroy,
            &storage,
            &Block<double>::save,
            &Block<double>::load);
    diy::ContiguousAssigner   assigner(world.size(), -1);   // number of blocks set by read_blocks()

     // read MFA model
    diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<double>::load);
    int nblocks = master.size();
    std::cout << nblocks << " blocks read from file "<< infile << "\n";



 // initialize DIY
    diy::FileStorage storage2("./DIY.YYYYYY"); // used for blocks to be moved out of core
    diy::Master      master2(world,
            -1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy,
            &storage2,
            &Block<real_t>::save,
            &Block<real_t>::load);
    diy::ContiguousAssigner   assigner2(world.size(), -1);   // number of blocks set by read_blocks()

    int nblocks2 = 0;

    if(func_type==2)
    {
        diy::io::read_blocks(infile2.c_str(), world, assigner2, master2, &Block<real_t>::load);
        nblocks2 = master2.size();
        std::cout << nblocks2 << " blocks read from file "<< infile2 << "\n";
    }



    master.foreach([&](Block<double>* b, const diy::Master::ProxyWithLink& cp)
    {
        
        Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;

        std::vector<VectorX<double>> point;
        std::vector<double> error;

        read_bin(root_file,point);

        switch (func_type)
        {
        case 0:
            getError(b,point,error,point_func_value);
            break;
        
        case 1:
            get_rv_error(b,point,error);
            break;

        case 2:
        {
             master2.foreach([&](Block<real_t>* b2, const diy::Master::ProxyWithLink& cp2)
            {
                get_js_error(b,b2,point,error);

            });
        }
            break;  
        }

   
        

        double min_error = *std::min_element(error.begin(), error.end());
        double max_error = *std::max_element(error.begin(), error.end());
        double sum_error = std::accumulate(error.begin(), error.end(), 0.0);
        double average_error = sum_error / error.size();

        printf("min_error: %e, max_error: %e, average_error: %e\n", min_error, max_error, average_error);

        std::ofstream file(func_error_file.c_str(), std::ios::binary);
        if (file.is_open()) {
            file.write(reinterpret_cast<const char*>(error.data()), error.size() * sizeof(double));
            file.close();
        } else {
            std::cout << "Failed to open file for writing." << std::endl;
        }
       

            

        // }

    });

    return 0;

}