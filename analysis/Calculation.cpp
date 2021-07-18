#include "Calculation.hpp"
Calculation::Calculation(InputPack& input)
{
    input.params().readString("name", KeyType::Required, name_);
    input.addCalculation(name_, this);
    return;
}
Calc_Nv::Calc_Nv(InputPack& input):Calculation{input}
{
    input.params().readString("probe_volume", KeyType::Required, pv_name_);
    auto pv_pointer = input.findProbeVolume(pv_name_);

}
void Calc_Nv::calculate(const Box& box){
  float sum = 0.0;
  for(int i = 0; i < box.atoms.size(); i++){
    sum += pv_->compute(box.atoms[i].x);
  }
  return;
}