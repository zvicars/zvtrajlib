#pragma once
#include "../ProbeVolume.hpp"
class PV_DiscreteRect : public ProbeVolume{
public: 
  PV_DiscreteRect(InputPack& input);
  virtual double compute(Vec3<double> position);
private:
  Vec<double> axis_ranges_;
};