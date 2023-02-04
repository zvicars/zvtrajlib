#include "PV_Simple.hpp"
PV_DiscreteRect::PV_DiscreteRect(std::array<double,2> xrange, std::array<double,2> yrange, std::array<double,2> zrange)
{
  for(int i = 0; i < 2; i++){
    axis_ranges_[i*3] = xrange[i];
    axis_ranges_[i*3+1] = yrange[i];
    axis_ranges_[i*3+2] = zrange[i];
  }
  return;
}

double PV_DiscreteRect::compute(Vec3<double> position){
  for(int i = 0; i < 3; i++){
    if( position[i] < axis_ranges_[i] || position[i] > axis_ranges_[i+3] ) return 0.0; 
  }
  return 1.0;
}