#include "Calc_Nv_wFit.hpp"

Calc_Nv_wFit::Calc_Nv_wFit(InputPack& input) : Calc_Nv{input} {
  std::string calcname;
  input.params().readString("calculation", KeyType::Required, calcname);
  auto calc_pointer = input.findCalculation(calcname);
  FANCY_ASSERT(calc_pointer != 0, "Failed to find specified calculation.");  
  calc_ = dynamic_cast<Calc_1D_Density*>(calc_pointer);
  if(calc_ == nullptr){
    calc2_ = dynamic_cast<Calc_1D_DensityFull*>(calc_pointer);
    FANCY_ASSERT(calc2_ != nullptr, "Dynamic cast failed in Calc_Nv_wFit");
    mode_ = 1;
  }
  else mode_ = 0;
  dir_ =  0;
  input.params().readFlag("direction", KeyType::Optional, dir_);
  return;
}

void Calc_Nv_wFit::calculate(){
  if(!doCalculate()) return;
  int dim;
  double x;
  if(mode_ == 0){
    dim = calc_->get_dim();
    x = calc_->get_x();
  }
  else{
    dim = calc2_->get_dim();
    x = calc2_->get_x();   
  }
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

