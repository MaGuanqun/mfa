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


#include "EigenMatrixIO.h"
#include"jacobi_set_h_critical_point.h"
#include "save_control_data.hpp"
#include "critical_point_jacobi_set.hpp"

#include"../critical_point/find_unique_root.h"


#include <chrono>


#include <tbb/tbb.h>
#include <atomic>

#include "../critical_point/index_of_critical_point.h"

using namespace std;

int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file
    string inControlPoint = "derivative_control_point.dat";
    string outfile = "derivative_root_point.dat";


    string infile2 = "approx2.mfa";               // diy input file
    string inControlPoint2 = "derivative_control_point2.dat";

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help

    int initial_point_level=1;
    // get command line arguments
    opts::Options ops;
    string input_shrink_ratio = "0-1-0-1";

    double same_root_epsilon = SAME_ROOT_EPSILON;
    double root_finding_epsilon = ROOT_FINDING_EPSILON;
    int max_itr=70;
    double function_value_threshold = 1e-6;

    //ops >> opts::Option('d', "deriv",   deriv,   " which derivative to take (1 = 1st, 2 = 2nd, ...)");
    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('g', "infile2",  infile2,  " diy input file name");
    ops >> opts::Option('i', "inControlPoint",  inControlPoint,  " diy input derivative control point file name");\
    ops >> opts::Option('j', "inControlPoint2",  inControlPoint2,  " diy input derivative control point file name");
    // ops >> opts::Option('a', "partial", partial, " dimension of 1 partial derivative only");
    ops >> opts::Option('l', "initial_point_level", initial_point_level, " 2^n+1 number of initial points");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('o', "outfile", outfile,    " control points of first derivative");
    ops >> opts::Option('e', "same_root_epsilon", same_root_epsilon, " epsilon for same root");
    ops >> opts::Option('s', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-z1-z2-...\"");
    ops >> opts::Option('t', "root_finding_epsilon", root_finding_epsilon, " root_finding_epsilon");
    ops >> opts::Option('m', "max_itr", max_itr, " max_itr");
    ops >> opts::Option('v', "function_value_threshold", function_value_threshold, " function_value_threshold");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // echo args
    fprintf(stderr, "\n--------- Input arguments ----------\n");
    // cerr <<
    //     "deriv = "    << deriv << endl;
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

    std::vector<std::vector<std::vector<std::vector<double>>>> geo_control_point; //[deriv][vars][dom][...]
    std::vector<std::vector<MatrixX<double>>> sci_deriv_control_points;//[vars][partial_deriv][...]
    save_control_points::load_control_points(inControlPoint.c_str(),geo_control_point,sci_deriv_control_points);



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


    // std::cout<<sci_deriv_control_points[0][0]<<std::endl;
    // std::cout<<"=="<<std::endl;
    // std::cout<<sci_deriv_control_points[0][1]<<std::endl;


    // read MFA model
    diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
    int nblocks = master.size();
    std::cout << nblocks << " blocks read from file "<< infile << "\n";



    //compute knots span [ui, ui+1) that could contain a critical point
    double decode_time = MPI_Wtime();


    std::vector<std::vector<VectorX<double>>> root(master.size());//[blocks,...]
    std::vector<std::vector<VectorX<double>>> domain_root(master.size()); //[blocks,]
    std::vector<std::vector<VectorX<double>>> domain_root_value_0(master.size()); //[blocks,]       
    //std::vector<std::vector<int>> multiplicity_root(master.size());
    int index = 0;

    std::vector<std::vector<int>> root_duplicated_number(master.size());

    std::vector<std::vector<VectorX<double>>> value(master.size());//[blocks,...]
    std::vector<std::vector<VectorX<double>>>deriv_value(master.size());//[blocks,...]
    
    std::vector<std::vector<int>> critical_point_type(master.size());

    int root_num=0;

    atomic_size_t filtered_out_num;
    filtered_out_num.store(0);

    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        {


        master2.foreach([&](Block<real_t>* b2, const diy::Master::ProxyWithLink& cp2)
        {


                Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;
                std::vector<std::vector<VectorXi>> selected_span;//[vars,span index]
                auto start_time = std::chrono::high_resolution_clock::now();


                double min_ = local_domain_range.minCoeff();
                auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
                VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;
                same_root_epsilon *= min_;
                VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());
                double min_span_size = Span_size.minCoeff()/10.0;
                if(min_span_size<same_root_epsilon)
                {
                    same_root_epsilon = min_span_size;
                }

                selected_span.clear();

                if(2*b->core_maxs.size()!=shrink_ratio.size())
                {
                    for(int i=0;i<b->core_maxs.size() - shrink_ratio.size()/2;i++)
                    {
                        shrink_ratio.push_back(0.0);
                        shrink_ratio.push_back(1.0);
                    }
                }

                for(int i=0;i<b->core_maxs.size();++i)
                {
                    shrink_ratio[2*i]= round((b->input->ndom_pts(i)-1)*shrink_ratio[2*i])/(b->input->ndom_pts(i)-1);
                    shrink_ratio[2*i+1]=round((b->input->ndom_pts(i)-1)*shrink_ratio[2*i+1])/(b->input->ndom_pts(i)-1);
                }    

                valid_span_jacobi_set::compute_valid_span(sci_deriv_control_points,b,selected_span,shrink_ratio);

                // std::cout<<"test span num "<<selected_span[0].size()<<std::endl;
                // for(int i=0;i<selected_span[0].size();++i)
                // {
                //     std::cout<<selected_span[0][i].transpose()<<std::endl;
                // }

                auto cpt_extract_start_time = std::chrono::high_resolution_clock::now();
                std::vector<int>multi_root_span; 
                // Eigen::VectorXd weights=Eigen::VectorXd::Ones(b->mfa->var(0).tmesh.tensor_prods[0].ctrl_pts.rows());

                std::vector<std::vector<VectorXd>> root_temp_record(selected_span[0].size());
                std::vector<int>original_root_size(selected_span[0].size(),0);

                tbb::affinity_partitioner ap;

                atomic_ulong itr_num;


                tbb::parallel_for(tbb::blocked_range<size_t>(0,selected_span[0].size()),
                [&](const tbb::blocked_range<size_t>& range)
                {
                    size_t filter_out=0;
                    for(auto i=range.begin();i!=range.end();++i)
                    {
                        size_t block_itr_num=0;
                        // std::vector<VectorXd> root_span;
                        jacobi_set_h_critical_point::newtonSolve(b,selected_span,root_temp_record[i],i,original_root_size[i],filter_out,same_root_epsilon,block_itr_num,root_finding_epsilon,max_itr,b2);
                        itr_num.fetch_add(block_itr_num);
                        // if(!root_span.empty())
                        // {
                        //     root_ori.insert(root_ori.end(),root_span.begin(),root_span.end());

                        //     // std::cout<<"checkl "<<root_span.size()<<" "<<root_span[0]<<std::endl;
                        // }
                        if(i%1000==0)
                        {
                            std::cout<<"finish find span "<<i<<std::endl;
                        }
                        // break;
                    }

                    filtered_out_num.fetch_add(filter_out);
                },ap               
                );


                std::cout<<"finish find all roots"<<std::endl;


                 for(auto i=0;i<selected_span[0].size();++i)
                 {
                    if(!root_temp_record[i].empty())
                    {
                        root[index].insert(root[index].end(),root_temp_record[i].begin(),root_temp_record[i].end());
                    }
                 }


                for(auto m=original_root_size.begin();m<original_root_size.end();++m)
                {
                    root_num+=(*m);
                }

                // std::cout<<"root size== "<<root_ori.size()<<std::endl;
                // std::cout<<root_ori[0]<<std::endl;
                // for(int i=0;i<root_ori.size();i++)
                // {
                //     std::cout<<"test "<<root_ori[i].transpose()<<std::endl;
                // }
                auto cpt_extract_end_time = std::chrono::high_resolution_clock::now();
                spatial_hashing::find_all_unique_root(root[index], domain_root[index],same_root_epsilon);                
                std::cout<<"root num before removing depulication "<<root_num<<std::endl;
                std::cout<<"root size "<<root[index].size()<<" "<<domain_root[index].size()<<std::endl;

                // std::cout<<"root size "<<root_ori.size()<<" "<<root[index].size()<<std::endl;
                // jacobi_set_h_critical_point::convertFromDomain(domain_root[index],root[index],b->core_mins,local_domain_range);

                
                root[index].clear();

                jacobi_set_h_critical_point::filterRoot(b,domain_root[index],root[index],domain_root_value_0[index],function_value_threshold, b2);



                // jacobi_set_h_critical_point::get_critical_point_index(b->mfa, b->vars[0].mfa_data,root_value_0[index],critical_point_type[index],HESSIAN_DET_EPSILON);

                // std::cout<<"critical point number for Jacobi set "<<root_value_0[index].size()<<std::endl;

                // jacobi_set_h_critical_point::getDerivative(b->mfa, b->vars[0].mfa_data,root[index],deriv_value[index]);
                // double result_gradient_average = jacobi_set_h_critical_point::getAverageNorm(deriv_value[index]);

                // std::cout<<"average gradient norm "<<result_gradient_average<<std::endl;

                auto end_time = std::chrono::high_resolution_clock::now();

                // std::vector<int> critical_point_index_;
                // critical_point_index::critical_point_index(b,root[index],critical_point_index_,HESSIAN_DET_EPSILON);
                // critical_point_index::output_index_summary(critical_point_index_);

                // std::vector<VectorXd> test;
                // jacobi_set_h_critical_point::find_all_unique_root(domain_root[index], test,same_root_epsilon,root_duplicated_number[index]);

                // std::cout<<"double check root number "<<test.size()<<std::endl;

                auto duration1=std::chrono::duration_cast<std::chrono::microseconds>(cpt_extract_end_time - cpt_extract_start_time);
                std::cout<<"cpt extract running time, millisecond : "<<duration1.count()/1000<<std::endl;
                
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                std::cout<<"root finding running time, millisecond : "<<duration.count()/1000<<std::endl;
                //get function value

                int initial_point_num = (b->mfa->var(0).p+VectorXi::Ones(b->mfa->var(0).p.size())).prod();
                double average = (double)(itr_num.load())/(double)(initial_point_num*selected_span[0].size());
                std::cout<<"Newton's method average itr num "<<average<<std::endl;

                index++;            
            });  
                
        });



    
    //save roots to a file
    std::vector<MatrixXd> root_matrix(domain_root_value_0.size());
    int domain_size, value_size;
    for(int i=0;i<domain_root_value_0.size();i++)
    {
        domain_size = root[i][0].size();
        value_size = root[i][0].size();
        root_matrix[i].resize(domain_root_value_0[i].size(),domain_size+value_size);
        for(int j=0;j<domain_root_value_0[i].size();j++)
        {
            root_matrix[i].block(j,0,1,domain_size) = root[i][j].transpose();
            root_matrix[i].block(j,domain_size,1,value_size) = domain_root_value_0[i][j].transpose();
            // root_matrix[i](j,domain_size+value_size) = 0;//(double)critical_point_type[i][j];
        }
    }

    // std::cout<<root_matrix[0]<<std::endl;

    
    std::cout<<"critical point number: "<< root_matrix[0].rows()<<std::endl;
    // std::cout<<root_matrix[0]<<std::endl;     
    // std::cout<<root_matrix[0]<<std::endl;

    utility::writeMatrixVector(outfile.c_str(),root_matrix);

    std::cout<<"filtered out number "<<filtered_out_num<<std::endl;

    // std::vector<MatrixXd> deplicate(domain_root.size());
    // for(int i=0;i<domain_root.size();i++)
    // {
    //     deplicate[i].resize(root_duplicated_number[i].size(),1);

    //     std::cout<<"duplicate root size "<<root_duplicated_number[i].size()<<std::endl;

    //     for(int j=0;j<root_duplicated_number[i].size();j++)
    //     {
    //         deplicate[i](j,0) = root_duplicated_number[i][j];
    //     }
    // }

    // string duplicate_outfile = "degree_of_"+outfile;
    // utility::writeMatrixVector(duplicate_outfile.c_str(),deplicate);


    // for(int i=0;i<deriv_value[0].size();i++)
    // {
    //     for(int j=0;j<deriv_value[0][i].size();++j)
    //     {
    //         std::cout<<deriv_value[0][i][j]<<std::endl;
    //     }
    // }

    // std::cout<<deriv_value[0][1].size()<<std::endl;


    // decode_time = MPI_Wtime() - decode_time;

    // std::cout<<selected_span.size()<<std::endl;

    // std::cout<<"root "<<std::endl;
    // for(auto i=0;i<root[0].size();++i)
    // {
    //     std::cout<<domain_root[0][i]<<" "<<multiplicity_root[0][i]<<std::endl;
    // }

    // Eigen::MatrixXd root_point(domain_root[0].size(),2);
    // root_point.setZero();
    // for(int i=0;i<domain_root[0].size();++i)
    // {
    //     root_point(i,0)=domain_root[0][i];
    // }
    // EigenMatrixIO::write_binary(outfile.c_str(),root_point);


    // std::vector<int>control_index;
    // for(auto i=selected_span.begin();i<selected_span.end();++i)
    // {
    //     for(auto j=0;j<4;++j)
    //     {
    //         control_index.emplace_back(*i-j);
    //     }
    // }

    // std::string na="selected_control_point.dat";
    // Eigen::MatrixXd used_control_point(control_index.size(),2);
    // for(auto i=0; i<control_index.size();++i)
    // {
    //     used_control_point.row(i)=deriv_control_point.row(control_index[i]);
    // }

    // EigenMatrixIO::write_binary(na.c_str(),used_control_point);

}

