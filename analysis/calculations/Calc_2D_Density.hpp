#pragma once
#include "Calculation.hpp"
#include <cmath>
struct GridDataExportPack{
  int integration_direction;
  std::vector<std::vector<double> > grid_data;
  std::array<int,2> grid_size, original_grid_size, grid_offset;
  std::array<double,2> real_offset, grid_spacing;
};

class Calc_2D_Density : public Calculation{
public:
  Calc_2D_Density(InputPack& input);
  virtual void calculate();
  virtual void update();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  GridDataExportPack getGridSubsection(std::array<int, 2> xmin, std::array<int,2> xmax);
private:
  void add_gaussian(double x_in, double y_in);
  std::array<int, 2> get_axes(){
    if(aligned_axis_ == 0){
      return {1, 2};
    }
    if(aligned_axis_ == 1){
      return {0, 2};
    } 
    if(aligned_axis_ == 2){
      return {0, 1};
    }
    throw 1;
    return {0,0};
  }
  int getBin(double x, int dimension){
    int bin;
    double dim_size = box_size_[dimension];
    if(x >= dim_size)  x-=dim_size;
    if(x < 0)  x+=dim_size;
    bin = floor(x/grid_spacing_[dimension]);
    return bin;
  }
  void putInBin(Vec3<double> pos){
    int idx = gridIndex(getBin(pos[ axes_[0] ], 0), getBin(pos[ axes_[1] ], 1) );
    grid_density_[idx] += 1.0;
    return;
  }
  int gridIndex(int x, int y){
    return npoints_[1]*x + y;
  }
  Vec<double> grid_density_;
  std::array<int,2> npoints_, axes_;
  std::array<double,2> grid_spacing_, box_size_, average_grid_spacing_;
  int aligned_axis_, frame_counter_; //x, y, or z as 0, 1, or 2, determines which dimension will be averaged out
  //this one depends on 3 atomgroups of equal size
  AtomGroup* atom_group_;
  bool initialized_, coarseGrain;
  double sigma_;
};