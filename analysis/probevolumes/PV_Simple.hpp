#pragma once
#include "ProbeVolume.hpp"
class PV_DiscreteRect : public ProbeVolume{
public: 
  PV_DiscreteRect(InputPack& input);
  PV_DiscreteRect(std::array<double,2> xrange, std::array<double,2> yrange, std::array<double,2> zrange);
  virtual double compute(Vec3<double> position);
private:
  Vec<double> axis_ranges_;
};