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


#include <tbb/tbb.h>
#include <atomic>

#include "find_isocontour.h"
#include"span_filter.h"
#include "find_root.h"
#include "../critical_point/find_all_root.h"
#include "find_jacobi_set.h"
#include "find_root_jacobi_set.h"

#include "connect_rv_graph.h"
#include "../critical_point/find_unique_root.h"

// #include "transfer_data.h"
// #include <fstream>
// #include <iostream>
// #include <Eigen/Dense>


using namespace std;

struct root_info
{
    std::vector<VectorX<real_t>> roots;
    std::vector<real_t> function_value;
    size_t valid_span_index;
};

// namespace {
//     tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
// }


int main(int argc, char** argv)
{

    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file
    string infile2 = "approx2.mfa";               // diy input file
    // string inControlPoint = "derivative_control_point.dat";
    string inCriticalPoint = "critical_point.dat";
    string inCriticalPoint2 = "critical_point2.dat";
    string root_file = "root_point.dat";

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;
    string input_point_string = "0-1-0-1";
    int max_step = 5000; //max point number in one isocontour

    double shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";

    string isosurface_file = "isocontour.dat";
    string ridge_valley_graph_file = "ridge_valley_graph.dat";
    double connect_traces_threshold = 2.0;

    string valley_file = "valley.dat";
    string ridge_file = "ridge.dat";

    string pseudo_valley_file = "pseudo_valley.dat";
    string pseudo_ridge_file = "pseudo_ridge.dat";


    double tolerance = 1e-4; // same_root_epsilon
    double step_size = 5.0;
    double epsilon = 1e-3; // decide if the result is close enough to the initail point and stop iteration

    double trace_split_grad_square_threshold = 1e-4; //check if need to split the trace

    int correction_max_itr = 30;

    // double center_threshold_square = epsilon*epsilon;
    double trace_threshold_square = step_size*step_size; //check duplication
    double threshold_connect_traces_square = 1e-4;

    int seperate_ridge_valley = 0;

    double point_func_value = 0.0; //value of the point need to be extracted.
    // double same_root_epsilon = SAME_ROOT_EPSILON;

    real_t root_finding_epsilon = 1e-8;

    real_t threshold_distance_stop_tracing = 0.1;

    string input_shrink_ratio = "0-1-0-1";

    int set_value_zero = 1;

    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('g', "infile2",  infile2,  " diy input file name");
    ops >> opts::Option('j', "inCriticalPoint",  inCriticalPoint,  " diy input critical point file name");
    ops >> opts::Option('c', "inCriticalPoint2",  inCriticalPoint2,  " diy input critical point file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('r', "root_file", root_file,    " root");
    ops >> opts::Option('t', "isosurface_file", isosurface_file, " file name of isosurface points");
    ops >> opts::Option('b', "ridge_valley_graph_file", ridge_valley_graph_file, " file name of ridge valley graph");
    ops >> opts::Option('e', "same_root_epsilon", epsilon, " epsilon for reaching start point");
    ops >> opts::Option('i', "input_point",    input_point_string,       " input point, by \"x1-x2-y1-y2-z1-z2-...\"");
    ops >> opts::Option('z', "step_size",    step_size,       " step size");
    ops >> opts::Option('v', "point_func_value",    point_func_value,       " value of the point need to be extracted");
    ops >> opts::Option('a', "connect_traces_threshold",    connect_traces_threshold,       " use adaptive RKF45 or not");
    ops >> opts::Option('m', "seperate_ridge_valley",    seperate_ridge_valley,       " seperate ridge valley in isocontour");
    ops >> opts::Option('p', "trace_split_grad_square_threshold",    trace_split_grad_square_threshold,       "determine split a trace");

    ops >> opts::Option('d', "ridge_file",    ridge_file,       "ridge point file");
    ops >> opts::Option('l', "valley_file",    valley_file,       "valley point file");

    ops >> opts::Option('s', "pseudo_ridge_file",    pseudo_ridge_file,       "pseudo ridge point file");
    ops >> opts::Option('q', "pseudo_valley_file",    pseudo_valley_file,       "pseudo valley point file");

    ops >> opts::Option('x', "root_finding_epsilon",    root_finding_epsilon,       "first root finding epsilon");

    ops >> opts::Option('y', "threshold_distance_stop_tracing",    threshold_distance_stop_tracing,       "when the distance from previous point is smaller than this*step_size, stop tracing");
    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");
    ops >> opts::Option('o', "set_value_zero",    set_value_zero,       " set the value of the point to zero");



    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }


    std::istringstream iss(input_shrink_ratio);
    std::vector<double> shrink_ratio;
    double number;
    std::string token;
    while (std::getline(iss, token, '-')) {
        std::istringstream tokenStream(token);
        if (tokenStream >> number) {
            shrink_ratio.push_back(number);
        }
    }  

    trace_threshold_square = step_size*step_size;
    std::cout<<"trace_threshold_square "<<trace_threshold_square<<" "<<step_size<<std::endl;

    epsilon = step_size;


    threshold_connect_traces_square = 1.21*step_size*step_size;

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

     // read MFA model
    diy::io::read_blocks(infile2.c_str(), world, assigner2, master2, &Block<real_t>::load);
    int nblocks2 = master2.size();
    std::cout << nblocks2 << " blocks read from file "<< infile2 << "\n";



    // std::vector<double> input_point_;
    // double number;
    // std::string token;
    // std::istringstream iss(input_point_string);
    // while (std::getline(iss, token, '-')) {
    //     std::istringstream tokenStream(token);
    //     if (tokenStream >> number) {
    //         input_point_.push_back(number);
    //     }
    // }  
    // VectorXd input_point(input_point_.size());
    // for(int i=0;i<input_point_.size();i++)
    // {
    //     input_point(i)=input_point_[i];
    // }
    // std::vector<int> sample_point_number_;
    // int number2;
    // std::istringstream iss2(input_sample_point_number);
    // while (std::getline(iss2, token, '-')) {
    //     std::istringstream tokenStream(token);
    //     if (tokenStream >> number2) {
    //         sample_point_number_.push_back(number2);
    //     }
    // } 
    // VectorXi sample_point_number(sample_point_number_.size());
    // for(int i=0;i<sample_point_number_.size();i++)
    // {
    //     sample_point_number(i)=sample_point_number_[i];
    // }

    // std::vector<real_t> values;
    // values.push_back(point_func_value);


    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {
        
        // Eigen::VectorXd weights=Eigen::VectorXd::Ones(b->vars[0].mfa_data->tmesh.tensor_prods[0].ctrl_pts.rows());
        Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;
        double min_ = local_domain_range.minCoeff();
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        // std::cout<<"span_num "<<span_num.transpose()<<std::endl;
        
        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());

        // std::cout<<"span size "<<Span_size.transpose()<<std::endl;
        double min_span_size = Span_size.minCoeff()/step_size;
        step_size = min_span_size;

        
        tolerance = 0.5*step_size; // same_root_epsilon
        epsilon = step_size; // decide if the result is close enough to the initail point and stop iteration
        trace_threshold_square = step_size*step_size; //check duplication
        // std::cout<<"step "<<step_size<<std::endl;
        // std::cout<<"trace_threshold_square__"<<trace_threshold_square<<std::endl;
        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);

        threshold_connect_traces_square = connect_traces_threshold*connect_traces_threshold*step_size*step_size;
        master2.foreach([&](Block<real_t>* b2, const diy::Master::ProxyWithLink& cp2)
        {
            std::vector<std::vector<VectorX<real_t>>> root; //the inner vector store the root in a span
            std::vector<std::vector<real_t>> function_value;
            std::vector<VectorX<real_t>> root_record;
            std::vector<real_t> record_function_value;
            std::vector<size_t> valid_span_index;
            std::vector<TraceInSpan<double>> traces_in_span;     

            int spanned_block_num =span_num.prod();

            std::vector<TraceInSpan<double>> ridge_valley_in_span;

            // std::cout<<"start finding root"<<std::endl;
            


            tbb::enumerable_thread_specific<std::vector<root_info>> local_root;

            // std::vector<std::vector<VectorX<real_t>>> root_block; //the inner vector store the root in a span
            // std::vector<std::vector<real_t>> function_value_block;
            // std::vector<size_t> valid_span_index_block;

            auto start_time = std::chrono::high_resolution_clock::now();

            tbb::affinity_partitioner ap;

            tbb::parallel_for(tbb::blocked_range<size_t>(0,spanned_block_num),
            [&](const tbb::blocked_range<size_t>& range)
            {

                auto& root_thread = local_root.local();

                for(auto i=range.begin();i!=range.end();++i)
                {
                    std::vector<VectorX<real_t>> root_block;
                    std::vector<real_t> function_value_block;
                    root_block.reserve(16);
                    function_value_block.reserve(16);
                    if(find_root_jacobi_set::root_finding(b,b2,root_block,number_in_every_domain,i, root_finding_epsilon,function_value_block,tolerance))
                    {
                        root_info span_root;
                        span_root.roots = std::move(root_block);
                        span_root.function_value = std::move(function_value_block);
                        span_root.valid_span_index = i;

                        root_thread.emplace_back(std::move(span_root));

                    }
                }
            },ap

            );

            size_t total_size = 0;
            for (const auto& thread_indices : local_root)
                total_size += thread_indices.size();

            root.reserve(total_size);
            function_value.reserve(total_size);
            valid_span_index.reserve(total_size);

            for (auto& thread_root_block : local_root)
            {
                for(auto& root_info : thread_root_block)
                {
                    root.emplace_back(std::move(root_info.roots));
                    valid_span_index.emplace_back(root_info.valid_span_index);
                    function_value.emplace_back(std::move(root_info.function_value));
                }
            }
                            


            auto root_end_time = std::chrono::high_resolution_clock::now();

            // std::cout<<"finish finding root"<<root.size()<<std::endl;

            int root_num=0;
            for(auto i=0;i<root.size();++i)
            {
                root_num+=root[i].size();
            }

            // std::cout<<"total root number "<< root_num<<std::endl;

            ridge_valley_in_span.resize(valid_span_index.size());

            jacobi_set::find_jacobi_set(step_size,max_step,b,b2,root, valid_span_index,trace_threshold_square,ridge_valley_in_span,root_finding_epsilon,trace_split_grad_square_threshold,correction_max_itr,threshold_distance_stop_tracing*threshold_distance_stop_tracing);

            // std::cout<<"finish tracing"<<std::endl;

            // std::cout<<"total points "<<num<<std::endl;

            jacobi_set::getFunctionValue(b,ridge_valley_in_span);
            // std::cout<<"getting function "<<std::endl;



            std::vector<VectorXd> ori_critical_points;
            connect_rv_graph::read_critical_point(inCriticalPoint,ori_critical_points,true);
            connect_rv_graph::read_critical_point(inCriticalPoint2,ori_critical_points,true);



            std::vector<VectorXd> critical_points;
            spatial_hashing::find_all_unique_root(ori_critical_points, critical_points,step_size);



            connect_rv_graph::connect_ridge_valley_graph(ridge_valley_in_span,critical_points,b,spanned_block_num, span_num,
            valid_span_index,
            threshold_connect_traces_square, ridge_valley_graph_file,step_size*step_size,set_value_zero);
            

            
            auto rv_graph_end_time = std::chrono::high_resolution_clock::now();


            // std::vector<MatrixXd> record_rv_surface(1);
            // find_isocontour::write_result_to_matrix(ridge_valley_in_span,record_rv_surface[0],b->core_mins,local_domain_range);

            // std::cout<<record_rv_surface[0].rows()<<" "<<record_rv_surface[0].cols()<<std::endl;

            // utility::writeMatrixVector(ridge_valley_graph_file.c_str(),record_rv_surface);


            std::cout<<"root finding running time, millisecond : "<< std::chrono::duration_cast<std::chrono::microseconds>(root_end_time - start_time).count()/1000<<std::endl;
            std::cout<<"ridge valley graph running time, millisecond : "<< std::chrono::duration_cast<std::chrono::microseconds>(rv_graph_end_time - root_end_time).count()/1000<<std::endl;
            std::cout<<"total running time, millisecond : "<< std::chrono::duration_cast<std::chrono::microseconds>(rv_graph_end_time - start_time).count()/1000<<std::endl;
            

            // jacobi_set::test_points_feature(b->mfa,b->vars[0].mfa_data,b2->mfa,b2->vars[0].mfa_data,ridge_valley_in_span,b->core_mins,local_domain_range);

            jacobi_set::get_js_error(b,b2,ridge_valley_in_span);

        });
            

            

        // }

    });


}