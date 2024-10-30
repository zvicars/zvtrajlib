#pragma once
#include "ProbeVolume.hpp"
class PV_INDUSRect : public ProbeVolume{
public: 
  PV_INDUSRect(InputPack& input);
  PV_INDUSRect(std::array<double,2> xrange, std::array<double,2> yrange, std::array<double,2> zrange);
  virtual double compute(Vec3<double> position);
private:
  std::array<double,6> axis_ranges_;
  double sigma_, cutoff_;
};