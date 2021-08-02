#include "VoxelGrid.hpp"
inline double heaviside(double x){
	if(x <= 0) return 0;
	return 1.0; 
}
inline double phi(double x, double sigma, double prefactor){
	return prefactor*(exp(-x*x/(2*sigma*sigma)))*heaviside(2.0*sigma-fabs(x));
}


inline double h_x(double x, double xmin, double xmax, double sigma, double xc){
    double eval = 0.0;
    double k, k1, k2, invk;
    double sigma2 = sigma*sigma;
    k = sqrt(2*M_PI)*sigma*erf(xc/(sqrt(2)*sigma)) - 2*xc*exp(-(xc*xc)/(2*sigma2));
    invk = 1/k;
    k1 = invk*sqrt(M_PI/2)*sigma;
    k2 = invk*exp(-(xc*xc)/(2*sigma2));
    eval = (k1 * erf((xmax-x)/(sqrt(2)*sigma)) - k2*(xmax-x) - 0.5)*heaviside(xc - fabs(xmax-x)) + heaviside(xc+xmax-x);
    return eval;
}


VoxelGrid::VoxelGrid(Vec3<int> size, Vec3<double> box_size, double density, double sigma, double isovalue) //must be orthorhombic box
{
  initialize(size, box_size, density, sigma, isovalue);
  return;
}
void VoxelGrid::initialize(Vec3<int> size, Vec3<double> box_size, double density, double sigma, double isovalue){
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
    grid_density_.resize(dim_x,std::vector<std::vector<double> >(dim_y,std::vector<double>(dim_z,0.0)));
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
        double xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = ix * grid_spacing_[0]
        xmax = xmin + grid_spacing_[0];
        ymin = iy * grid_spacing_[1];
        ymax = ymin + grid_spacing_[1];
        zmin = iz * grid_spacing_[2];
        zmax = zmin + grid_spacing_[2];
        grid_density_[idx][idy][idz] += (1.0/density_) * h_x(x, xmin, xmax, sigma_, 2.0*sigma_)*h_x(y, ymin, ymax, sigma_, 2.0*sigma_)*h_x(z, zmin, zmax, sigma_, 2.0*sigma_);
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