#include "VoxelGrid.hpp"
inline double heaviside(double x){
	if(x <= 0) return 0;
	return 1.0; 
}
inline double phi(double x, double sigma, double prefactor){
	return prefactor*(exp(-x*x/(2*sigma*sigma)))*heaviside(2.0*sigma-fabs(x));
}
VoxelGrid::VoxelGrid(Vec3<int> size, Vec3<double> box_size, double density, double sigma, double isovalue) //must be orthorhombic box
{
    int dx = size[0];
    int dy = size[1];
    int dz = size[2];
    resize_grid(dx, dy, dz);
    for(int i = 0; i < 3; i++){
        grid_spacing_[i] = box_size[i] / (double)(size[i]);
    }
    sigma_= sigma;
    density_ = density;
    isovalue_ = isovalue;
    prefactor_ = pow(2*M_PI*sigma_*sigma_, (-1.0/2.0));
    return;
}
int VoxelGrid::resize_grid(int dim_x, int dim_y, int dim_z)
{
    grid_density_.clear();
    std::vector<double> col(dim_z, 0.0);
    std::vector<std::vector<double> > plane(dim_y, col);
    grid_density_.resize(dim_x, plane);
    sz[0] = dim_x; sz[1]= dim_y; sz[2] = dim_z;
    return 1;
}
void VoxelGrid::add_gaussian(Vec3<double> x_in)
{
    int idx, idy, idz;
    double x = x_in[0];
    double y = x_in[1];
    double z = x_in[2];
    for(int ix = floor((x-2*sigma_)/grid_spacing_[0]); ix <= ceil((x+2*sigma_)/grid_spacing_[0]); ix++)
    for(int iy = floor((y-2*sigma_)/grid_spacing_[1]); iy <= ceil((y+2*sigma_)/grid_spacing_[1]); iy++)
    for(int iz = floor((z-2*sigma_)/grid_spacing_[2]); iz <= ceil((z+2*sigma_)/grid_spacing_[2]); iz++)
    {
        idx = ix; idy = iy; idz = iz;
        pbcidx(idx, idy, idz);
        grid_density_[idx][idy][idz] += (1.0/density_) * phi(ix*grid_spacing_[0] - x, sigma_, prefactor_)*phi(iy*grid_spacing_[1] - y, sigma_, prefactor_)*phi(iz*grid_spacing_[2]- z, sigma_, prefactor_);
    }
    return;
}
void VoxelGrid::pbcidx(int &x, int& y, int& z)
{
    if(x >= sz[0]) x = x%sz[0];
    else if(x < 0) x = sz[0] - (-x)%sz[0];
    if(y >= sz[1]) y = y%sz[1];
    else if(y < 0) y = sz[1] - (-y)%sz[1];
    if(z >= sz[2]) z = z%sz[2];
    else if(z < 0) z = sz[2] - (-z)%sz[2];
    return;
}