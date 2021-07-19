//Interface to several marching cubes algorithms
//Core functionality is to take a voxel grid and convert it to a series of triangles that can be analyzed/outputted
//Author: Zachariah Vicars
//Marching cubes algorithms come from various sources
#pragma once
#include <vector>
#include <array>
#include <string>
#include <cmath>

class VoxelGrid
{
public:
	template <class T>
	using Vec3 = std::array<T, 3>;
	VoxelGrid(Vec3<int> size, Vec3<double> box_size, double density, double sigma, double isovalue);
	int resize_grid(int dim_x, int dim_y, int dim_z, double ival = 0.0);
	Vec3<int> getSize() const{
		return sz;
	}
	Vec3<double> getLength() const{
		Vec3<double> ret;
		for(int i = 0; i < 3; i++){
			ret[i] = sz[i]*grid_spacing_;
		}
	}
	double getGridVal(int i, int j, int k) const{
		return grid_density_[i][j][k];
	}
	double getIsovalue() const{
		return isovalue_;
	}
	double get_gs() const{
		return grid_spacing_;
	}
    void pbcidx(int &x, int& y, int& z);
    void add_gaussian(double x, double y, double z, int idx1);

private:
	std::vector<std::vector<std::vector<double > > > grid_density_;
	Vec3<int> sz; 
	double grid_spacing_;
	double prefactor_, sigma_, isovalue_;
};