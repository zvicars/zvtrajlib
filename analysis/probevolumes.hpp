#pragma once
#include "../interface/datatypes.hpp"
#include "parser.hpp" //to be able to register itself

class ProbeVolume{ //probe volumes allow you to specify a particle position and receive a floating point number corresponding to whether the particle is in the volume or not
public:
  ProbeVolume(std::string input);
  virtual double compute(Vec3<double> position) = 0;
private:
  std::string name_;
};

class PV_DiscreteRect : public ProbeVolume{
public: 
  PV_DiscreteRect(std:string input):ProbeVolume{input};
  double calculate(Vec3<double> position);
private:
  Vec3<double> lower_corner_;
  Vec3<double> upper_corner_;
};

