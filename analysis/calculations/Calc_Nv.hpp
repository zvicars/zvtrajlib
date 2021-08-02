#pragma once
#include "Calculation.hpp"
#include <cmath>
#include "../helper/make_histogram.hpp"
class Calc_Nv : public Calculation{
public:
  Calc_Nv(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:
  Vec<double> time_vec_;
  Vec<int> step_vec_;
  Vec<double> count_vec_;
  double value_, mean_, var_;
  std::string pv_name_;
  ProbeVolume* pv_;
  AtomGroup* atom_group_; //these are the atoms for which calculations will be performed
};