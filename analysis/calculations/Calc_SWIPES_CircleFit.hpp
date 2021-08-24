#pragma once
#include "Calculation.hpp"
#include "../helper/make_histogram.hpp"
#include "Calc_2D_Density.hpp"
#include <cmath>

//fits a circle to a 2d grid density, accepts a left-most boundary 
class Calc_SWIPES_CircleFit : public Calculation{
public:
  Calc_SWIPES_CircleFit(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:
  ProbeVolume* pv_; //will only consider atoms within the bounded region
  Calc_2D_Density* calc_;
  //atomgroup to be considered
  AtomGroup* atom_group_;
  double density_;
  //subsection of the density grid to consider
  Vec3<int> xmin_, xmax_;
  //tells the algorithm which direction to iterate over to place vertices, either 0 or 1, for the respective coordinate of the Calc_2D_Density grid
  int normal_direction_;
  //placed vertices
  Vec<Vec3<double> > vertex_positions_;
  //circlefit parameters, x0, y0, R
  Vec3<double> params_;
  //whether or not the parameter will be participating in the fit
  Vec3<bool> fix_params_;
};