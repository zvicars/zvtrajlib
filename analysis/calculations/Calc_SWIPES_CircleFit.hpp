#pragma once
#include "Calculation.hpp"
#include "Calc_2D_Density.hpp"
#include <cmath>

//fits a circle to a 2d grid density, accepts a left-most boundary 
class Calc_SWIPES_CircleFit : public Calculation{
public:
  Calc_SWIPES_CircleFit(InputPack& input);
  virtual void finalOutput();
protected:
  Calc_2D_Density* calc_;
  ProbeVolume* pv_; //will only consider atoms within the bounded region
  //atomgroup to be considered
  AtomGroup* atom_group_;
  double density_;
  //subsection of the density grid to consider
  std::array<int,2> xmin_, xmax_;
  //tells the algorithm which direction to iterate over to place vertices, either 0 or 1, for the respective coordinate of the Calc_2D_Density grid
  int normal_direction_, parallel_direction_;
  //placed vertices
  Vec<std::array<double,2> > vertex_positions_;
  //circlefit parameters, x0, y0, R
  Vec3<double> params_;
  //whether or not the parameter will be participating in the fit
  Vec3<bool> fix_params_;
};