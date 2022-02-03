#pragma once
#include "Calc_Histogram.hpp"
#include <cmath>
class Calc_Nv : public Calculation_Histogram{
public:
  Calc_Nv(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void update(){
    Calculation::update();
    return;
  }
protected:
  Vec<double> time_vec_;
  Vec<int> step_vec_;
  Vec<double> count_vec_;
  double value_, mean_, var_;
  std::string pv_name_;
  ProbeVolume* pv_;
  AtomGroup* atom_group_; //these are the atoms for which calculations will be performed
  int dump_frame_ = -1;
};