#include "MarchingCubesInterface.hpp"
#include <sstream>
namespace golosio{
#include "./golosio/TriangulateGlobal.hpp"
}
void printSTL(const Mesh& mesh, std::string& frame)
{
	std::stringstream ofile;
	ofile << "solid " << "frame" << "\n";
	for(int i = 0; i < mesh.ntri; i++)
	{
        int vidx1 = mesh.triangles[i].indices[0];
        int vidx2 = mesh.triangles[i].indices[1];
        int vidx3 = mesh.triangles[i].indices[2];
        //just using the first vertex normal as the face normal, lazy but shouldn't affect much
        ofile << "facet normal " << mesh.normals[vidx1][0] << "  " << mesh.normals[vidx1][1] << "  " << mesh.normals[vidx1][2] << "\n"; 
		ofile << "    outer loop\n";
		ofile << "vertex " << mesh.vertices[vidx1][0] << " " << mesh.vertices[vidx1][1] << " " << mesh.vertices[vidx1][2] << "\n";
 		ofile << "vertex " << mesh.vertices[vidx2][0] << " " << mesh.vertices[vidx2][1] << " " << mesh.vertices[vidx2][2] << "\n";
		ofile << "vertex " << mesh.vertices[vidx3][0] << " " << mesh.vertices[vidx3][1] << " " << mesh.vertices[vidx3][2] << "\n";       
		ofile << "    endloop\n";
		ofile << "endfacet\n";
	}
	ofile << "endsolid " << frame << std::endl;
	frame = ofile.str();
  return;
}

void marchingCubes(std::string type, const VoxelGrid& v, Mesh& output_mesh)
{
	int nx = v.getSize()[0];
	int ny = v.getSize()[1];
	int nz = v.getSize()[2];
	double box_vec[3] = {v.getLength()[0], v.getLength()[1], v.getLength()[2]};
	std::vector<float> volume_data(nx*ny*nz, 0);
	int iterator = 0;
	for(int k = 0; k < nz; k++)
	for(int j = 0; j < ny; j++)
	for(int i = 0; i < nx; i++)
	{
		volume_data[iterator] = v.getGridVal(i, j, k);
		iterator++;
	}
	float* vertex_data; 
	int* triangle_data;
	int nvtx; 
	int ntri; 
	int n[3] = {nx, ny, nz};

	float s[3] = {(float)v.get_gs()[0], (float)v.get_gs()[1], (float)v.get_gs()[2]}; 
	float thresh = (float)v.getIsovalue(); 
	golosio::TriangulateGlobal g1;
	g1.Triangulate(&volume_data[0], &vertex_data, &triangle_data, n, s, &nvtx, &ntri, thresh, 0);
	std::vector<Vec3<double> > vertices_(nvtx);
	std::vector<Vec3<double> > gradients_(nvtx);
	std::vector<Triangle> triangles_(ntri);
	for(int i = 0; i < nvtx; i++){
			for(int j = 0; j < 3; j++){
					vertices_[i][j] = vertex_data[6*i + j] + 0.5*v.getLength()[j];
					gradients_[i][j] = vertex_data[6*i + j + 3];   
			}
	}
	for(int i = 0; i < ntri; i++){
        Vec3<int> triangle_indices_;
        for(int j = 0; j < 3; j++){
            triangle_indices_[j] = triangle_data[3*i+j];   
        }
        triangles_[i].indices = triangle_indices_;
	}
    output_mesh.normals = gradients_;
    output_mesh.vertices = vertices_;
    output_mesh.triangles = triangles_;
    output_mesh.nvtx = nvtx;
    output_mesh.ntri = ntri;

	free(vertex_data);
	free(triangle_data);
	return;
}
