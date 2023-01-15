#include "Calc_Nv_wFit_IP.hpp"
#include "../../tools/pbcfunctions.hpp"
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
  double com_dx = calc_->get_com_dx();
  float sum = 0.0;
  for(int i = 0; i < atom_group_->getIndices().size(); i++){
    int idx = atom_group_->getIndices()[i];
    auto pos = box->atoms[idx].x;
    pos[dim] = wrapNumber(pos[dim] + com_dx, box->boxvec[dim][dim]);
    if((pos[dim] >= x1 && pos[dim] <= x2) != dir_){
      sum += pv_->compute(pos);
    }
  }
  value_ = sum;
  if(doTimeseries || doHistogram){
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  return;
}

