//Interface to several marching cubes algorithms
//Core functionality is to take a voxel grid and convert it to a series of triangles that can be analyzed/outputted
//Author: Zachariah Vicars
//Marching cubes algorithms come from various sources
#pragma once
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
class VoxelGrid
{
public:
	template <class T>
	using Vec3 = std::array<T, 3>;
	VoxelGrid(){
		return;
	}
	VoxelGrid(Vec3<int> size, Vec3<double> box_size, double density, double sigma, double isovalue);
	void initialize(Vec3<int> size, Vec3<double> box_size, double density, double sigma, double isovalue);
	int resize_grid(int dim_x, int dim_y, int dim_z);
	Vec3<int> getSize() const{
		return sz;
	}
	Vec3<double> getLength() const{
		Vec3<double> ret;
		for(int i = 0; i < 3; i++){
			ret[i] = (sz[i])*grid_spacing_[i];
		}
		return ret;
	}
	void setLength(Vec3<double> length){
		for(int i = 0; i < 3; i++){
			grid_spacing_[i] = length[i]/(double)(sz[i]);
		}
		return;
	}
	double getGridVal(int i, int j, int k) const{
		return grid_density_[i][j][k];
	}
	double getIsovalue() const{
		return isovalue_;
	}
	Vec3<double> get_gs() const{
		return grid_spacing_;
	}
	void sumInPlace(VoxelGrid* vg){
		for(int i = 0; i < sz[0]; i++){
			for(int j = 0; j < sz[1]; j++){
				for(int k = 0; k < sz[2]; k++)
				{
					std::cout << i << j << k << std::endl;
					grid_density_[i][j][k] += vg->getGridVal(i,j,k);
				}
			}
		}
		return;
	}
	void scalarMult(double rhs) // compound assignment (does not need to be a member,
	{                           // but often is, to modify the private members)
		for(int i = 0; i < sz[0]; i++){
			for(int j = 0; j < sz[1]; j++){
				for(int k = 0; k < sz[2]; k++){
					grid_density_[i][j][k] *= rhs;
				}
			}
		}
		return;
	}
	void pbcidx(int &x, int& y, int& z);
	void add_gaussian(Vec3<double> x_in);
	void clear(){
		for(int i = 0; i < sz[0]; i++){
			for(int j = 0; j < sz[1]; j++){
				for(int k = 0; k < sz[2]; k++)
				{
					grid_density_[i][j][k] = 0.0; 
				}
			}
		}
		return;
	}
private:
	std::vector<std::vector<std::vector<double > > > grid_density_;
	Vec3<int> sz; 
	Vec3<double> grid_spacing_;
	double prefactor_, sigma_, isovalue_, density_;
};