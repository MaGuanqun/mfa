#include <vector>
#include <iostream>
#include <Eigen/Dense>

namespace write_to_ply
{
    void writeHeader(std::ofstream& file, size_t numVertices, size_t numFaces) {
        file << "ply\n";
        file << "format binary_little_endian 1.0\n"; // Use binary_little_endian or binary_big_endian as needed
        file << "element vertex " << numVertices << "\n";
        file << "property float x\n";
        file << "property float y\n";
        file << "property float z\n";
        file << "element face " << numFaces << "\n";
        file << "property list uchar int vertex_index\n"; // uchar: number of vertices per face, int: vertex indices
        file << "end_header\n";
    }

    template<typename T>
    void writeVertices(std::ofstream& file, const std::vector<std::vector<T>>& vertex_domain,
    std::vector<T>& pt_data) {
    size_t i=0;

        float index_1_f;
        float index_0_f;
        float value;
        for (auto index_1:vertex_domain[1]) {
            index_1_f = index_1;
            for(auto index_0:vertex_domain[0]) {
                index_0_f = index_0;
                value = pt_data[i];
                file.write(reinterpret_cast<const char*>(&index_0_f), sizeof(index_0_f));
                file.write(reinterpret_cast<const char*>(&index_1_f), sizeof(index_1_f));
                file.write(reinterpret_cast<const char*>(&value), sizeof(value));
                i++;
            }
        }
    }

    void writeFaces(std::ofstream& file, VectorXi ndom_pts) {
        unsigned char numVertices = 3;
        for (unsigned int i = 0; i < ndom_pts(1) - 1; ++i) {
            for (unsigned int j = 0; j < ndom_pts(0) - 1; ++j) {
                unsigned int v0 = i * ndom_pts(0) + j;
                unsigned int v1 = v0 + 1;
                unsigned int v2 = v0 + ndom_pts(0);
                unsigned int v3 = v2 + 1;
                file.write(reinterpret_cast<const char*>(&numVertices), sizeof(numVertices));
                file.write(reinterpret_cast<const char*>(&v0), sizeof(unsigned int));
                file.write(reinterpret_cast<const char*>(&v1), sizeof(unsigned int));
                file.write(reinterpret_cast<const char*>(&v2), sizeof(unsigned int));
                file.write(reinterpret_cast<const char*>(&numVertices), sizeof(numVertices));
                file.write(reinterpret_cast<const char*>(&v1), sizeof(unsigned int));
                file.write(reinterpret_cast<const char*>(&v3), sizeof(unsigned int));
                file.write(reinterpret_cast<const char*>(&v2), sizeof(unsigned int));
            }
        }
    }

    template<typename T>
    void write_ply(std::string& filename, const std::vector<std::vector<T>>& vertex_domain,
    std::vector<T>& pt_data)
    {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }
        int num_triangle =2* (vertex_domain[0].size()-1) * (vertex_domain[1].size()-1);
        writeHeader(file, vertex_domain[0].size()*vertex_domain[1].size(), num_triangle);

        writeVertices(file, vertex_domain, pt_data);
        VectorXi ndom_pts(2);
        ndom_pts << vertex_domain[0].size(), vertex_domain[1].size();
        writeFaces(file, ndom_pts);
        file.close();
    }
}