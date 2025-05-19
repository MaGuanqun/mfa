#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "EigenMatrixIO.h"
#include"find_unique_root.h"
#include    <diy/master.hpp>
#include    <diy/io/block.hpp>

#include    "opts.h"

#include    "block.hpp"

#include    "mfa/mfa.hpp"

#include "utility_function.h"

void csv_processing(std::string& infile, std::vector<Eigen::VectorXd>& points)
{
    std::string line, cell;

    // Open the CSV file
    std::ifstream file(infile);

    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        exit(1);
    }

    // Skip the first line (title line)
    if (!getline(file, line)) {
        std::cerr << "Error reading the title line from the file." << std::endl;
        exit(1);
    }



    double x_max = -1e10;
    double x_min = 1e10;
    double y_max = -1e10;
    double y_min = 1e10;

    // Read the file line by line
    while (getline(file, line)) {
        std::stringstream lineStream(line);
        std::string type, is_boundary;
        double x, y, z;
        // Read the type and is_boundary columns
        getline(lineStream, type, ',');
        getline(lineStream, is_boundary, ',');

        // Skip the boundary points

        if (is_boundary == "1") {
            continue;
        }
        // Read the x, y, z coordinates
        getline(lineStream, cell, ',');
        x = std::stod(cell);
        getline(lineStream, cell, ',');
        y = std::stod(cell);
        getline(lineStream, cell);
        z = std::stod(cell);

        Eigen::VectorXd point(2);
        point << x, y;
        // Add the point to the vector
        points.push_back(point);
    }

    file.close();
}

void find_min_max(int argc, char** argv, string& input_mfa_file,VectorXd& min, VectorXd& max)
{
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD


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



    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {  
        min = b->core_mins;
        max = b->core_maxs;
         std::cout<<" mfa min max are: "<<b->core_mins.transpose()<<" "<<b->core_maxs.transpose()<<std::endl;
    });
}




void proceed_find_unique_root(std::vector<Eigen::VectorXd>& points, VectorXd& min, VectorXd& max,std::vector<Eigen::VectorXd>& unique_points, double same_root_epsilon)
{
       // normalize point


    // for (auto& point : points) {
        
    //     point(0) = (point(0) - min(0)) / (max(0) - min(0));
    //     point(1) = (point(1) - min(1)) / (max(1) - min(1));
    // }

    // for (auto& point : points) {
    //     printf("point(0) %f point(1) %f\n",point(0),point(1));
    // }

    // Find the unique points
   

    spatial_hashing::find_all_unique_root(points,unique_points,same_root_epsilon);

    // Write the unique points to a file
    // std::cout<<"unique_points.size() "<<unique_points.size()<<std::endl;
    // std::cout<<"points.size() "<<points.size()<<std::endl;

    // std::sort(unique_points.begin(), unique_points.end(), compareVectors);
    // for (auto& point : unique_points) {
    //     std::cout<<point.transpose()<<std::endl;
    // }
    // for (auto& point : unique_points) {
    //     printf("point(0) %f point(1) %f\n",point(0),point(1));
    // }

    // for (auto& point : unique_points) {
    //     point(0) = point(0) * (max(0) - min(0))+min(0);
    //     point(1) = point(1) * (max(1) - min(1))+min(1);
    // }

    std::cout<<"points.size() "<<points.size()<<std::endl;
    std::cout<<"unique_points.size() "<<unique_points.size()<<std::endl;
    // printf("min(0) %f min(1) %f\n",min(0),min(1));
    // printf("max(0) %f max(1) %f\n",max(0),max(1));
    

}


void operate_based_on_extension(int argc, char** argv, string& input_mfa_file, std::string& filename,string& outfile, double value, double same_root_epsilon) {
    // Find the last dot in the filename
    size_t lastDotPos = filename.find_last_of(".");
    if (lastDotPos != std::string::npos) {
        // Extract the extension
        std::string extension = filename.substr(lastDotPos + 1);
        
        // // Convert the extension to lowercase for case-insensitive comparison (optional)
        // std::transform(extension.begin(), extension.end(), extension.begin(),
        //                [](unsigned char c) { return std::tolower(c); });

        // Perform operations based on the extension
        std::vector<Eigen::VectorXd> points;
        if (extension == "dat") {
      
        }
        else if (extension == "csv") {

            csv_processing(filename, points);

            
        }
        else {
            std::cout << "Unknown file extension: " << extension << std::endl;
        }

        VectorXd min, max;
        find_min_max(argc, argv, input_mfa_file, min, max);

        std::vector<Eigen::VectorXd> unique_points;
        proceed_find_unique_root(points,min,max,unique_points,same_root_epsilon);
        utility::saveToCSV(outfile, unique_points,value);
    }
}


int main(int argc, char** argv) {

    std::string infile = "name.csv";  

    opts::Options ops;

 
    bool help;  
    string input_mfa_file = "approx.mfa";
    string outfile = "root.csv";
    double value=0;
    double same_root_epsilon = SAME_ROOT_EPSILON;

    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('i', "input_mfa_file",      input_mfa_file,     " input mfa file name");
    ops >> opts::Option('o', "outfile", outfile,    " output csv file");
    ops >>opts::Option('v', "value",    value,    " function value");
    ops >>opts::Option('e', "same_root_epsilon",    same_root_epsilon,    " same root epsilon");

    if (!ops.parse(argc, argv) || help)
    {
        std::cout << ops;
        return 1;
    }

    #ifdef MFA_TBB
        std::cerr << "threading: TBB" << std::endl;
    #endif

    std::cerr << "infile = " << infile << std::endl;

    


    operate_based_on_extension(argc, argv, input_mfa_file, infile, outfile, value,same_root_epsilon);

 



}
