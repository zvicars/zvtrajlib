#include "Calc_Nv.hpp"
namespace CalculationRegistry {
static const Register<Calc_Nv>
  registerType("nv");
}
Calc_Nv::Calc_Nv(InputPack& input):Calculation{input}
{
  input.params().readString("probe_volume", KeyType::Required, pv_name_);
  auto pv_pointer = input.findProbeVolume(pv_name_);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer;

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");

  return;
}
void Calc_Nv::calculate(){
  if(!doCalculate()) return;
  float sum = 0.0;
  for(int i = 0; i < atom_group_->getIndices().size(); i++){
    int idx = atom_group_->getIndices()[i];
    sum += pv_->compute(box->atoms[idx].x);
  }
  value_ = sum;
  if(doOutput()) printOutput();
  return;
}
std::string Calc_Nv::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Nv" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}