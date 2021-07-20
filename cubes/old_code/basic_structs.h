#pragma once
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstdint>

inline std::string trim(const std::string &s)
{
   auto wsfront=std::find_if_not(s.begin(),s.end(),[](int c){return std::isspace(c);});
   auto wsback=std::find_if_not(s.rbegin(),s.rend(),[](int c){return std::isspace(c);}).base();
   return (wsback<=wsfront ? std::string() : std::string(wsfront,wsback));
}
struct atom
{
	double pos[3];
	std::string name;
};
double heaviside(double x)
{
	if(x <= 0) return 0;
	return 1.0; 
}

double phi(double x, double sigma, double prefactor)
{
	return 1.05*prefactor*(exp(-x*x/(2*sigma*sigma)))*heaviside(2.0*sigma-fabs(x));
}

class voxel_grid
{
private:
	std::vector<std::vector<std::vector<double > > > grid_density;
	int size_x, size_y, size_z;
	int frame_counter = 0;
	double grid_spacing;
	std::string atom_type; 
	std::vector<double> sigma_vec, density_vec, prefactor_vec;
public:
	voxel_grid();
	voxel_grid(int dx, int dy, int dz, double gs, std::vector<double> sigma1, std::vector<double> dens)
	{
		resize_grid(dx, dy, dz);
		size_x = dx;
		size_y = dy;
		size_z = dz;
		grid_spacing = gs;
		sigma_vec = sigma1;
		density_vec = dens;
		for(int i = 0; i < sigma1.size(); i++)
		{
			prefactor_vec.push_back(pow(2*M_PI*sigma_vec[i]*sigma_vec[i], (-1.0/2.0)));
		}
		return;
	}
	int resize_grid(int dim_x, int dim_y, int dim_z, double ival = 0.0)
	{
		grid_density.clear();
		std::vector<double> col(dim_z, ival);
		std::vector<std::vector<double> > plane(dim_y, col);
		grid_density.resize(dim_x, plane);
		size_x = dim_x; size_y = dim_y; size_z = dim_z;
		return 1;
	}
	double getlength(int dim)
	{
		if (dim == 0) return size_x * grid_spacing; 
		if (dim == 1) return size_y * grid_spacing;
		if (dim == 2) return size_z * grid_spacing;
		else return 0.0;
	}
	int getsize(int dim)
	{
		if (dim == 0) return size_x; 
		if (dim == 1) return size_y;
		if (dim == 2) return size_z;
		else return 0;		
	}
	void add_gaussian(double x, double y, double z, int idx1)
	{
		int idx, idy, idz;
		double sigma = sigma_vec[idx1];
		double prefactor = prefactor_vec[idx1];
		double density = density_vec[idx1];
		for(int ix = floor((x-2*sigma)/grid_spacing); ix <= ceil((x+2*sigma)/grid_spacing); ix++)
		for(int iy = floor((y-2*sigma)/grid_spacing); iy <= ceil((y+2*sigma)/grid_spacing); iy++)
		for(int iz = floor((z-2*sigma)/grid_spacing); iz <= ceil((z+2*sigma)/grid_spacing); iz++)
		{
			idx = ix; idy = iy; idz = iz;
			pbcidx(idx, idy, idz);
			grid_density[idx][idy][idz] += (1.0/density) * phi(ix*grid_spacing - x, sigma, prefactor)*phi(iy*grid_spacing - y, sigma, prefactor)*phi(iz*grid_spacing - z, sigma, prefactor);
		}
		return;
	}
	void pbcidx(int &x, int& y, int& z)
	{
		if(x >= size_x) x = x%size_x;
		else if(x < 0) x = size_x - (-x)%size_x;
		if(y >= size_y) y = y%size_y;
		else if(y < 0) y = size_y - (-y)%size_y;
		if(z >= size_z) z = z%size_z;
		else if(z < 0) z = size_z - (-z)%size_z;
		return;
	}
	double get_grid_val(int i, int j, int k)
	{
		return grid_density[i][j][k]/((double)frame_counter);
	}
	double get_gs()
	{
		return grid_spacing;
	}
	int it_counter()
	{
		frame_counter++;
		return frame_counter;
	}
};
