//save different data to a matrix
#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>


namespace transfer_data
{
    template<class T>
    void transfer(std::vector<std::vector<Eigen::VectorX<T>>>& isocontour_points_in_domain,std::vector<std::vector<Eigen::VectorXd>>& isocontour_func_value, Eigen::MatrixXd& record_iso_surface) {

        int row_num=0;
        std::vector<int> accumulate_number(isocontour_points_in_domain.size()+1,0); //record the start index of every vector if we combine them together
        for(int i=0;i <isocontour_points_in_domain.size();++i)
        {
            row_num+=isocontour_points_in_domain[i].size();
            accumulate_number[i+1]=accumulate_number[i]+isocontour_points_in_domain[i].size();
        }

        record_iso_surface.resize(row_num,isocontour_points_in_domain[0][0].size()+isocontour_func_value[0][0].size());
        for(int i = 0; i < isocontour_points_in_domain.size(); ++i)
        {
            int start_index = accumulate_number[i];
            for(int j=0;j<isocontour_points_in_domain[i].size();++j)
            {
                record_iso_surface.block(start_index+j,0,1,isocontour_points_in_domain[i][j].size())=isocontour_points_in_domain[i][j].transpose();
                record_iso_surface.block(start_index+j,isocontour_points_in_domain[i][j].size(),1,isocontour_func_value[i][j].size())=isocontour_func_value[i][j].transpose();
            }
        }

    }

    
}