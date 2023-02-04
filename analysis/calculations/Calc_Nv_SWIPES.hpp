#pragma once
#include "Calc_Histogram.hpp"
#include <cmath>
#include "../probevolumes/PV_Simple.hpp"
class Calc_Nv_SWIPES : public Calculation_Histogram{
public:
  Calc_Nv_SWIPES(InputPack& input);
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
  //this one will make its own probevolumes
  std::vector<PV_DiscreteRect> pv_set_;
  std::vector<std::vector<double> > ts_values_;
  std::vector<std::array<double, 2> > ranges_;
  AtomGroup* atom_group_; //these are the atoms for which calculations will be performed
};