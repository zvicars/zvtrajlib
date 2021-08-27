#include "Calc_Nv_wFit.hpp"

Calc_Nv_wFit::Calc_Nv_wFit(InputPack& input) : Calc_Nv{input} {
  std::string calcname;
  input.params().readString("calculation", KeyType::Required, calcname);
  auto calc_pointer = input.findCalculation(calcname);
  FANCY_ASSERT(calc_pointer != 0, "Failed to find specified calculation.");  
  calc_ = dynamic_cast<Calc_1D_Density*>(calc_pointer);
  dir_ =  0;
  input.params().readFlag("direction", KeyType::Optional, dir_);
  return;
}

void Calc_Nv_wFit::calculate(){
  if(!doCalculate()) return;
  int dim = calc_->get_dim();
  double x = calc_->get_x();
  float sum = 0.0;
  for(int i = 0; i < atom_group_->getIndices().size(); i++){
    int idx = atom_group_->getIndices()[i];
    if((box->atoms[idx].x[dim] <= x) != dir_) sum += pv_->compute(box->atoms[idx].x);
  }
  value_ = sum;
  if(doTimeseries || doHistogram){
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  return;
}

