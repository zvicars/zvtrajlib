#pragma once
#include "ProbeVolume.hpp"
class PV_Trivial : public ProbeVolume{
public: 
  PV_Trivial(InputPack& input);
  PV_Trivial();
  virtual double compute(Vec3<double> position);
private:
  std::array<double,6> axis_ranges_;
};