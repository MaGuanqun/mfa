#include<random>

#include"find_unique_root.h"
#include"find_all_root.h"
#include <chrono>

int main()
{

    const int numPoints = 100000;
    const double rangeMin = -5.0;
    const double rangeMax = 5.0;

    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(rangeMin, rangeMax);

    // Generate the points
    std::vector<Eigen::VectorX<double>> points;
    for (int n = 0; n < numPoints; ++n) {
        Eigen::VectorX<double> vec(3);
        vec<<dis(gen), dis(gen), dis(gen);
        points.push_back(vec);
    }

    // Introduce some duplications
    std::uniform_int_distribution<> indexDis(0, numPoints - 1);
    int numDuplications = 100;  // number of duplications

    Eigen::VectorX<double> turbu(3);
    turbu<<-6e-6, 0, 0;

    for (int i = 0; i < numDuplications; ++i) {
        int idx = indexDis(gen);
        points.push_back(points[idx]+turbu);
    }

    std::cout<<"original size "<<points.size()<<std::endl;

    std::vector<Eigen::VectorX<double>> unique_points_1;
    std::vector<int> duplicated_number;
    auto start_time1 = std::chrono::high_resolution_clock::now();
    find_all_roots::find_all_unique_root(points,unique_points_1,SAME_ROOT_EPSILON,duplicated_number);
    auto end_time1 = std::chrono::high_resolution_clock::now();
    std::vector<Eigen::VectorX<double>> unique_points_2;

    auto start_time2 = std::chrono::high_resolution_clock::now();
    spatial_hashing::find_all_unique_root(points,unique_points_2,SAME_ROOT_EPSILON);
    auto end_time2 = std::chrono::high_resolution_clock::now();

    std::cout<<"unique root kd tree "<<unique_points_1.size()<<" "<<unique_points_2.size()<<std::endl;

    std::cout<<"kd-tree time "<<std::chrono::duration_cast<std::chrono::microseconds>(end_time1 - start_time1).count()/1000.0<<std::endl;
    std::cout<<"spatial hashing time "<<std::chrono::duration_cast<std::chrono::microseconds>(end_time2 - start_time2).count()/1000.0<<std::endl;



}