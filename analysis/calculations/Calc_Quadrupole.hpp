#pragma once
#include "Calculation.hpp"
#include "Eigen/Eigen"
#include <cmath>
class Calc_Quadrupole : public Calculation{
public:
  Calc_Quadrupole(InputPack& input);
  virtual void calculate();
  virtual void update();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:
  AtomGroup* atom_group_;
  //expecting a sequential list of molecules, will modulo with molsize to get appropriate charge for atom
  int natoms_;
  Vec<double> masses_, charges_;
  Vec3<double> normal_; 
  Vec<double> 
};