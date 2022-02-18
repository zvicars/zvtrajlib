#include "ProbeVolume.hpp"

ProbeVolume::ProbeVolume(InputPack& input){
  input.params().readString("name", KeyType::Required, name_);
  input.addProbeVolume(name_, this);
  box = input.getBox();
  update_flag_ = 0;
  return;
}