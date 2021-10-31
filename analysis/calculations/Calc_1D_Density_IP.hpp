#pragma once
#include "Calculation.hpp"
#include <cmath>


class Calc_1D_Density_IP : public Calculation{
public:
  Calc_1D_Density_IP(InputPack& input);
  virtual void calculate();
  virtual void update();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  std::array<double,4> getFit(){
    if(!fitSigmoidal){
      std::cout << "1D density profile is not fitting a sigmoidal." << std::endl;
      throw 0;
    }
    return params_;
  }
  double get_x(){
    return params_[3];
  }
  double get_x2(){
    return params_[2];
  }
  int get_dim(){
    return dim_;
  }
private:
  int getBin(double x){
    int bin;
    double dim_size = box_size_;
    if(x >= dim_size)  x-=dim_size;
    if(x < 0)  x+=dim_size;
    bin = floor(x/grid_spacing_);
    return bin;
  }
  void putInBin(Vec3<double> pos){
    int idx = getBin(pos[ dim_ ]);
    grid_density_[idx] += 1.0;
    return;
  }
  void add_gaussian(double x_in);
  Vec<double> grid_density_, average_grid_density_;
  int dim_, npoints_, frame_counter_; //the axis that will be kept
  double grid_spacing_, box_size_, average_grid_spacing_;

  AtomGroup* atom_group_;

  //fitting options
  bool fitSigmoidal;
  std::array<int, 2> idx_range_;
  Vec<double> guess_; // calculated per-frame
  
  Vec<double> tvec_, frame_vec_;
  Vec< std::array<double,4> > fits_;
  std::array<double,4> params_;

  //coarse-graining
  bool coarseGrain;
  double sigma_;
};