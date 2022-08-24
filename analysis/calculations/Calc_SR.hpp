//Zachariah Vicars 7/7/2022
//Calc SR is an abstract class meant to calculate short-ranged interactions, the parent class
//will handle the core algorithms needed to compute short-ranged interactions while the child
//classes will override the "compute()" function to allow for different potential types

#pragma once
#include "Calc_Histogram.hpp"
#include "../../tools/cellgrid.hpp"
#include <unordered_map>
#include <cmath>
class Calc_SR : public Calculation_Histogram{
public:
  Calc_SR(InputPack& input);
  virtual void calculate();
  virtual double compute_potential(int idx1, int idx2) = 0;
  virtual void finalOutput();
  virtual void update();
protected:
  Vec<double> time_vec_;
  Vec<int> step_vec_;
  Vec<double> value_vec_;
  double value_, mean_, var_;
  AtomGroup* atom_group_; //these are the atoms for which calculations will be performed
  double cutoff_;
  Vec3<double> box_size_;
  //per atom value of potential, SUM(U_step) = 2.0*U_tot
  std::vector<double> U_step_;
};