#include "PV_Boolean.hpp"

PV_Boolean::PV_Boolean(InputPack& input):ProbeVolume{input}
{
  std::string pv1_name, pv2_name;
  input.params().readString("pv1", KeyType::Required, pv1_name);
  auto pv1_pointer = input.findProbeVolume(pv1_name);
  FANCY_ASSERT(pv1_pointer != 0, "Failed to find specified probe volume.");
  pv1_ = pv1_pointer;  
  input.params().readString("pv2", KeyType::Required, pv2_name);
  auto pv2_pointer = input.findProbeVolume(pv2_name);
  FANCY_ASSERT(pv2_pointer != 0, "Failed to find specified probe volume.");
  pv2_ = pv2_pointer;  
  std::vector<bool> logic;
  input.params().readVector("logic", KeyType::Required, logic);
  FANCY_ASSERT(logic.size() == 4, "Invalid logic gate specified. Need 4 entries defining the behavior of a logic gate with different inputs.");
  for(int i = 0; i < 4; i++){
    logic_chart_[i] = logic[i];
  }
  return;
}

double PV_Boolean::compute(Vec3<double> position){
  bool eval1 = pv1_->compute(position);
  bool eval2 = pv2_->compute(position);
  if(eval1){
    if(eval2){
      return logic_chart_[3];
    }
    else{
      return logic_chart_[2];
    }
  }
  else{
    if(eval2){
      return logic_chart_[1];
    }
    else{
      return logic_chart_[0];
    }
  }
  return 0;
}