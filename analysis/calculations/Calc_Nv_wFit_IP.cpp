#include "Calc_Nv_wFit_IP.hpp"

Calc_Nv_wFit_IP::Calc_Nv_wFit_IP(InputPack& input) : Calc_Nv{input} {
  std::string calcname;
  input.params().readString("calculation", KeyType::Required, calcname);
  auto calc_pointer = input.findCalculation(calcname);
  FANCY_ASSERT(calc_pointer != 0, "Failed to find specified calculation.");  
  calc_ = dynamic_cast<Calc_1D_Density_IP*>(calc_pointer);
  dir_ =  0;
  input.params().readFlag("direction", KeyType::Optional, dir_);
  return;
}

void Calc_Nv_wFit_IP::calculate(){
  if(!doCalculate()) return;
  if(!(calc_->hasCalculated())) calc_->calculate();
  int dim = calc_->get_dim();
  double x1 = calc_->get_x();
  double x2 = calc_->get_x2();
  float sum = 0.0;
  for(int i = 0; i < atom_group_->getIndices().size(); i++){
    int idx = atom_group_->getIndices()[i];
    if((box->atoms[idx].x[dim] >= x1 && box->atoms[idx].x[dim] <= x2) != dir_) sum += pv_->compute(box->atoms[idx].x);
  }
  value_ = sum;
  if(doTimeseries || doHistogram){
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  return;
}

