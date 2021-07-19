#pragma once
#include <array>
#include <vector>
#include "VoxelGrid.hpp"

template <class T>
using Vec3 = std::array<T, 3>;

struct Triangle{
    Vec3<int> indices; //contains either the real positions of the vertices or the indices, depending on whether you use an int or a double
    Vec3<Vec3<double > > normals; //contains the normal vector of a triangle
};
struct Mesh {
    std::vector<Vec3<double> > vertices;
    std::vector<Vec3<double> > normals;
    std::vector<Triangle> triangles;
    int nvtx, ntri;
};

//Main marching cubes function, will be able to switch between different algorithms
void marchingCubes(std::string type, const VoxelGrid& grid, Mesh& output_mesh);
void printSTL(const Mesh& mesh, std::string& frame);