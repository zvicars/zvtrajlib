#pragma once
#include "../Calculation.hpp"
#include "../helper/make_histogram.hpp"
#include <cmath>
class Calc_Angle : public Calculation{
public:
  Calc_Angle(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:
  Vec<Vec<double> > angle_distribution_;
  Vec<double> angle_vec_;
  Vec<int> step_vec_;
  Vec<double> time_vec_;
  double mean_, var_, value_;
  //this one depeends on 3 atomgroups of equal size
  Vec3<AtomGroup*> atom_groups_;
};