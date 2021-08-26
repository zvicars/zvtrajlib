#pragma once
#include "ProbeVolume.hpp"
#include "Eigen/Eigen"
class PV_Cylinder : public ProbeVolume{
public: 
  PV_Cylinder(InputPack& input);
  virtual double compute(Vec3<double> position);
private:
  Eigen::Vector3d position_;
  Eigen::Vector3d axis_;
  double radius_;
  double height_;
  Eigen::Matrix3d tmat_;
};