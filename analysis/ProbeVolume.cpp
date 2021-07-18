#include "ProbeVolume.hpp"
//place in registry

ProbeVolume::ProbeVolume(InputPack& input){
  input.params().readString("name", KeyType::Required, name_);
  input.addProbeVolume(name_, this);
  return;
}