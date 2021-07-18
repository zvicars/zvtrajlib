#pragma once
#include "ProbeVolume.hpp"

class PV_DiscreteRect : public ProbeVolume{
public: 
  PV_DiscreteRect(InputPack& input);
  double calculate(Vec3<double> position);
private:
  Vec3<double> axis_ranges_;
};
