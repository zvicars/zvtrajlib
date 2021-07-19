#include "VoxelGrid.hpp"
inline double heaviside(double x){
	if(x <= 0) return 0;
	return 1.0; 
}
inline double phi(double x, double sigma, double prefactor){
	return prefactor*(exp(-x*x/(2*sigma*sigma)))*heaviside(2.0*sigma-fabs(x));
}
VoxelGrid(Vec3<int> size, Vec3<double> box_size, double density, double sigma, double isovalue) //must be orthorhombic box
{
    resize_grid(dx, dy, dz);
    for(int i = 0; i < 3; i++){
        sz[i] = size[i];
        gs[i] = box_size[i] / (double)(size[i]-1);
    }
    grid_spacing_ = gs;
    sigma_= = sigma;
    density_ = density;
    isovalue_ = isovalue;
    prefactor_ = pow(2*M_PI*sigma_vec[i]*sigma_vec[i], (-1.0/2.0));
    return;
}
int VoxelGrid::resize_grid(int dim_x, int dim_y, int dim_z, double ival = 0.0)
{
    grid_density.clear();
    std::vector<double> col(dim_z, ival);
    std::vector<std::vector<double> > plane(dim_y, col);
    grid_density.resize(dim_x_, plane);
    sz[0] = dim_x; sz[1]= dim_y; sz[2] = dim_z;
    return 1;
}
void VoxelGrid::add_gaussian(double x, double y, double z)
{
    int idx, idy, idz;

    for(int ix = floor((x-2*sigma)/grid_spacing_); ix <= ceil((x+2*sigma)/grid_spacing_); ix++)
    for(int iy = floor((y-2*sigma)/grid_spacing_); iy <= ceil((y+2*sigma)/grid_spacing_); iy++)
    for(int iz = floor((z-2*sigma)/grid_spacing_); iz <= ceil((z+2*sigma)/grid_spacing_); iz++)
    {
        idx = ix; idy = iy; idz = iz;
        pbcidx(idx, idy, idz);
        grid_density[idx][idy][idz] += (1.0/density_) * phi(ix*grid_spacing_ - x, sigma_, prefactor_)*phi(iy*grid_spacing_ - y, sigma, prefactor_)*phi(iz*grid_spacing_- z, sigma_, prefactor_);
    }
    return;
}
void VoxelGrid::pbcidx(int &x, int& y, int& z)
{
    if(x >= size_x_) x = x%size_x_;
    else if(x < 0) x = size_x_ - (-x)%size_x_;
    if(y >= size_y_) y = y%size_y_;
    else if(y < 0) y = size_y_ - (-y)%size_y_;
    if(z >= size_z_) z = z%size_z_;
    else if(z < 0) z = size_z_ - (-z)%size_z_;
    return;
}