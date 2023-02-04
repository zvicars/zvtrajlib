#include "ProbeVolume.hpp"

ProbeVolume::ProbeVolume(InputPack& input){
  input.params().readString("name", KeyType::Required, name_);
  input.addProbeVolume(name_, this);
  box = input.getBox();
  update_flag_ = 0;
  return;
}

//default constructor for in-place construction of probe volume
ProbeVolume::ProbeVolume(std::string name){
  name_ = name;
  box = 0;
  update_flag_ = 0;
  return;
}
ProbeVolume::ProbeVolume(){
  name_ = "";
  box = 0;
  update_flag_ = 0;
  return;
}