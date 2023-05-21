#pragma once
#include "Calc_Histogram.hpp"
#include <cmath>
class Calc_HPolarity : public Calculation_Histogram{
public:
  Calc_HPolarity(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void update(){
    Calculation_Histogram::update();
    FANCY_ASSERT(h_group_->getIndexCount() == t_group_->getIndexCount(), "head and tail group indices are different sizes!");
    for(int i = 0; i < 3; i++){
      box_size_[i] = box->boxvec[i][i];
    }
    return;
  }
protected:
  Vec<double> time_vec_;
  Vec<int> step_vec_;
  Vec<double> count_vec_;
  double value_, mean_, var_;
  std::string pv_name_;
  ProbeVolume* pv_;
  AtomGroup* h_group_, * t_group_; //these are the atoms for which calculations will be performed
  std::array<double, 3> normal_, box_size_;
  int dump_frame_ = -1;
};