#include "Calc_Nv.hpp"
namespace CalculationRegistry {
static const Register<Calc_Nv>
  registerType("Nv");
}
Calc_Nv::Calc_Nv(InputPack& input):Calculation{input}
{
  input.params().readString("probe_volume", KeyType::Required, pv_name_);
  auto pv_pointer = input.findProbeVolume(pv_name_);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer; 
  return;
}
void Calc_Nv::calculate(const Box& box){
  float sum = 0.0;
  for(int i = 0; i < box.atoms.size(); i++){
    sum += pv_->compute(box.atoms[i].x);
  }
  value_ = sum;
  return;
}
std::string Calc_Nv::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Nv" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}