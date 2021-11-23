#pragma once
#include "ProbeVolume.hpp"
#include "../atomgroups/AtomGroup.hpp"
class PV_DynBox: public ProbeVolume{
public: 
  PV_DynBox(InputPack& input);
  virtual double compute(Vec3<double> position);
  virtual void update();
private:
  Vec<double> axis_ranges_;
  AtomGroup* atomgroup_; //needs atom positions to work
  Box* box_; //needs access to atom positions
  std::array<double, 3> xcom_; 
};