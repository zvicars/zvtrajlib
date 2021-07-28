#pragma once
#include "Calculation.hpp"
class Calc_Nv : public Calculation{
public:
  Calc_Nv(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
private:
  Vec<double> time;
  Vec<int> step;
  Vec<double> count;
  double value_;
  std::string pv_name_;
  ProbeVolume* pv_;
  AtomGroup* atom_group_; //these are the atoms for which calculations will be performed
};