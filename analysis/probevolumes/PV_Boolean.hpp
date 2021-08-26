#pragma once
#include "ProbeVolume.hpp"
class PV_Boolean : public ProbeVolume{
public: 
  PV_Boolean(InputPack& input);
  virtual double compute(Vec3<double> position);
private:
  ProbeVolume* pv1_;
  ProbeVolume* pv2_;
  //Defines what to do when..
  //[0] : A = 0, B = 0
  //[1] : A = 0, B = 1
  //[2] : A = 1, B = 0
  //[3] : A = 1, B = 1
  //AND 0 0 0 1
  //NOR 1 0 0 0
  //XOR 0 1 1 0
  //used because I need to be able to switch between AND and a non-commutative logic gate for SWIPES
  std::array<bool, 4> logic_chart_; 
};