#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include "EigenMatrixIO.h"
#include <filesystem>


namespace save_control_points
{
    template<class T>
    void save_control_points(const char* filename, std::vector<std::vector<std::vector<std::vector<T>>>>& geo_control_point, //[deriv][vars][dom][...]
            std::vector<std::vector<MatrixX<T>>>& sci_deriv_control_points) //[vars][partial_deriv][...]
    {


        // Create the directory path if it does not exist
        std::filesystem::path file_path(filename);
        auto parent_path = file_path.parent_path();
        if (!parent_path.empty() && !std::filesystem::exists(parent_path)) {
            std::filesystem::create_directories(parent_path);
        }
        //save vars, deriv, geo_domain, sci_domain
        //save geo_control, sci_control
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        int deriv, vars_num, domain_dimension;
        deriv= geo_control_point.size();//0,1st,2nd ...order derivative
        out.write((char*)(&deriv), sizeof(int));
        size_t size_in_domain;
        for(int i=0;i<deriv;i++)
        {
            vars_num = geo_control_point[i].size();
            out.write((char*)(&vars_num), sizeof(int));
            for(int j=0;j<vars_num;j++)
            {
                domain_dimension = geo_control_point[i][j].size();
                out.write((char*)(&domain_dimension), sizeof(int));
                for(int k=0;k<domain_dimension;k++)
                {
                    size_in_domain = geo_control_point[i][j][k].size();
                    out.write((char*)(&size_in_domain), sizeof(size_t));
                    out.write((char*)(geo_control_point[i][j][k].data()), size_in_domain*sizeof(T));
                }
            }
        }

        //save sci_control
        vars_num = sci_deriv_control_points.size();
        out.write((char*)(&vars_num), sizeof(int));
        int partial;
        size_t matrix_col, matrix_row;
        for(int i=0;i<vars_num;i++)
        {
            partial = sci_deriv_control_points[i].size();
            out.write((char*)(&partial), sizeof(int));
            for(int j=0;j<partial;j++)
            {
                matrix_row = sci_deriv_control_points[i][j].rows();
                matrix_col = sci_deriv_control_points[i][j].cols();
                out.write((char*)(&matrix_row), sizeof(size_t));
                out.write((char*)(&matrix_col), sizeof(size_t));
                out.write((char*)(sci_deriv_control_points[i][j].data()), matrix_row*matrix_col*sizeof(T));
            }
        }
		out.close();
    }

    template<class T>
    bool load_control_points(const char* filename, std::vector<std::vector<std::vector<std::vector<T>>>>& geo_control_point, //[deriv][vars][dom][...]
            std::vector<std::vector<MatrixX<T>>>& sci_deriv_control_points) //[vars][partial_deriv][...]
    {

        std::ifstream in(filename, std::ios::in | std::ios::binary);
		if (!in.good())
		{
			std::cout << "file not open" << std::endl;
			return false;
		}
        int deriv, vars_num, domain_dimension;
        in.read((char*)(&deriv), sizeof(int));
        geo_control_point.resize(deriv);
        size_t size_in_domain;
        for(int i=0;i<deriv;i++)
        {
            in.read((char*)(&vars_num), sizeof(int));
            geo_control_point[i].resize(vars_num);
            for(int j=0;j<vars_num;j++)
            {
                in.read((char*)(&domain_dimension), sizeof(int));
                geo_control_point[i][j].resize(domain_dimension);
                for(int k=0;k<domain_dimension;k++)
                {
                    in.read((char*)(&size_in_domain), sizeof(size_t));
                    geo_control_point[i][j][k].resize(size_in_domain);
                    in.read((char*)(geo_control_point[i][j][k].data()), size_in_domain*sizeof(T));
                }
            }
        }

        //read sci_control
        in.read((char*)(&vars_num), sizeof(int));
        sci_deriv_control_points.resize(vars_num);
        int partial;
        size_t matrix_col, matrix_row;
        for(int i=0;i<vars_num;i++)
        {
            in.read((char*)(&partial), sizeof(int));
            sci_deriv_control_points[i].resize(partial);
            for(int j=0;j<partial;j++)
            {
                in.read((char*)(&matrix_row), sizeof(size_t));
                in.read((char*)(&matrix_col), sizeof(size_t));
                sci_deriv_control_points[i][j].resize(matrix_row,matrix_col);
                in.read((char*)(sci_deriv_control_points[i][j].data()), matrix_row*matrix_col*sizeof(T));

            }
        }

		in.close();
		return true;
    }

}