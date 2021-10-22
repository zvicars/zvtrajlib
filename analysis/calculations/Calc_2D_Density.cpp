#include "Calc_2D_Density.hpp"

/*
class Calc_2D_Density : public Calculation{
public:
  Calc_2D_Density(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:

  Vec<Vec<double> > grid_density_;
  std::array<int,2> npoints_;
  std::array<double,2> grid_spacing_;
  int aligned_axis_; //x, y, or z as 0, 1, or 2, determines which dimension will be averaged out
  //this one depends on 3 atomgroups of equal size
  AtomGroup* atom_group_;
*/

inline double heaviside(double x){
	if(x <= 0) return 0;
	return 1.0; 
}

inline double h_x(double x, double xmin, double xmax, double sigma, double xc){
    double eval = 0.0;
    double k, k1, k2, invk;
    double sigma2 = sigma*sigma;
    k = sqrt(2*M_PI)*sigma*erf(xc/(sqrt(2)*sigma)) - 2*xc*exp(-(xc*xc)/(2*sigma2));
    invk = 1/k;
    k1 = invk*sqrt(0.5*M_PI*sigma2);
    k2 = invk*exp(-0.5*(xc*xc)/sigma2);
    eval = (k1 * erf((xmax-x)/(sqrt(2)*sigma)) - k2*(xmax-x) - 0.5)*heaviside(xc - fabs(xmax-x))
    + (k1 * erf((x-xmin)/(sqrt(2)*sigma)) - k2*(x-xmin) - 0.5)*heaviside(xc - fabs(x-xmin))
    + heaviside(xc + 0.5*(xmax-xmin) - fabs(x - 0.5*(xmin+xmax)));
    return eval;
}

Calc_2D_Density::Calc_2D_Density(InputPack& input) : Calculation{input} {
  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  std::vector<int> points;
  input.params().readVector("npoints", KeyType::Required, points);
  FANCY_ASSERT(points.size() == 2, "Invalid number of elements provided in " + name_ + ", npoints. Expected 2.");
  for(int i = 0; i < 2; i++){
    npoints_[i] = points[i];
  }
  input.params().readNumber("axis", KeyType::Required, aligned_axis_);
  FANCY_ASSERT(aligned_axis_ >= 0 && aligned_axis_ < 3, "Invalid dimension provided in " + name_ + ". Expected 0 for x, 1 for y, or 2 for z.");
  
  coarseGrain = 0;
  input.params().readFlag("smear", KeyType::Optional, coarseGrain);
  if(coarseGrain){
    input.params().readNumber("sigma", KeyType::Optional, sigma_);
  }

  axes_ = get_axes();
  frame_counter_ = 0;
  initialized_ = 0;

}
void Calc_2D_Density::update(){
  Calculation::update();
  Vec3<double> box_size;
  for(int i = 0; i < 3; i++){
    box_size[i] = box->boxvec[i][i];
  }
  if(!initialized_){
    grid_density_.resize(npoints_[0]*npoints_[1], 0.0);
    average_grid_spacing_ = {0};
    initialized_ = 1;
  }
  for(int i = 0; i < 2; i++){
    grid_spacing_[i] = box_size[axes_[i]] / (double)(npoints_[i]);
    box_size_[i] = box_size[axes_[i]];
  }
  return;
}
void Calc_2D_Density::calculate(){
  if(!doCalculate()) return;
  for(int i = 0; i < atom_group_->getIndices().size(); i++ ){
    int idx = atom_group_->getIndices()[i];
    auto position = box->atoms[idx].x;
    if(coarseGrain) add_gaussian(position[axes_[0]], position[axes_[1]]);
    else putInBin(position);
  }
  average_grid_spacing_[0] += grid_spacing_[0];
  average_grid_spacing_[1] += grid_spacing_[1];
  frame_counter_++;
}

void Calc_2D_Density::add_gaussian(double x_in, double y_in)
{
    double x = x_in;
    double y = y_in;
    int lxmin = floor((x-2*sigma_)/grid_spacing_[0]);
    int lxmax = ceil((x+2*sigma_)/grid_spacing_[0]); 
    int lymin = floor((y-2*sigma_)/grid_spacing_[1]);
    int lymax = ceil((y+2*sigma_)/grid_spacing_[1]);    

    #pragma omp parallel for collapse(2)
    for(int ix = lxmin; ix <= lxmax; ix++)
    {
      for(int iy = lymin; iy <= lymax; iy++)
      {
        int idx, idy;
        idx = ix; idy = iy;
        if(idx >= grid_density_.size()) idx -= npoints_[0];
        else if(idx < 0) idx += npoints_[0];
        if(idy >= grid_density_.size()) idy -= npoints_[1];
        else if(idy < 0) idy += npoints_[1];       
        double xmin, xmax, ymin, ymax;
        xmin = ix * grid_spacing_[0];
        xmax = xmin + grid_spacing_[0];
        ymin = iy * grid_spacing_[1];
        ymax = ymin + grid_spacing_[1];        
        grid_density_[gridIndex(getBin(x,0), getBin(x,1))] += h_x(x, xmin, xmax, sigma_, 2.0*sigma_)*h_x(y, ymin, ymax, sigma_, 2.0*sigma_);
      }
    }
    return;
}

std::string Calc_2D_Density::printConsoleReport(){
  return "";
}
void Calc_2D_Density::finalOutput(){
  average_grid_spacing_[0] *= 1.0/(double)frame_counter_;
  average_grid_spacing_[1] *= 1.0/(double)frame_counter_;
  double invterm =  1.0 / ((double)frame_counter_ * average_grid_spacing_[0] * average_grid_spacing_[1] * box->boxvec[aligned_axis_][aligned_axis_]);
  for(int i = 0; i < grid_density_.size(); i++){
    grid_density_[i] *= invterm;
  }
  std::ofstream ofile(base_ + "_2D_Density.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_)
  
  for(int i = 0; i < npoints_[0]; i++){
    for(int j = 0; j < npoints_[1]; j++){
      ofile << i*grid_spacing_[0] << "     " << j*grid_spacing_[1] << "     " << grid_density_[gridIndex(i, j)] << "\n";
    }
  }
  return;
}

GridDataExportPack Calc_2D_Density::getGridSubsection(std::array<int, 2> xmin, std::array<int,2> xmax){
  FANCY_ASSERT(xmin[0] > 0 && xmin[0] <= xmax[0], "Invalid xmin data provided to getGridSubsection in Calc_2D_Density.");
  FANCY_ASSERT(xmin[1] > 0 && xmin[1] <= xmax[1], "Invalid xmin data provided to getGridSubsection in Calc_2D_Density.");
  FANCY_ASSERT(xmax[0] < npoints_[0], "Invalid xmax data provided to getGridSubsection in Calc_2D_Density.");
  FANCY_ASSERT(xmax[1] < npoints_[1], "Invalid xmax data provided to getGridSubsection in Calc_2D_Density.");

  GridDataExportPack g1;
  g1.original_grid_size = npoints_;
  g1.grid_size = {xmax[0] - xmin[0] + 1, xmax[1] - xmin[1] + 1};
  g1.grid_offset = xmin;
  g1.real_offset = {xmin[0]*grid_spacing_[0], xmin[1]*grid_spacing_[1]};
  g1.grid_data.resize(xmax[0] - xmin[0] + 1, std::vector<double>(xmax[1] - xmin[1] + 1, 0.0));
  g1.grid_spacing = grid_spacing_;
  g1.integration_direction = aligned_axis_;
  for(int i = xmin[0]; i <= xmax[0]; i++){
    for(int j = xmin[1]; j <= xmax[1]; j++){
      int di = i - xmin[0];
      int dj = j - xmin[1];
      g1.grid_data[di][dj] = grid_density_[gridIndex(i,j)];
    }
  }
  return g1;
}