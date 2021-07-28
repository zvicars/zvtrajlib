#include "factory.hpp"
ProbeVolume* ProbeVolume_Factory(std::string key, InputPack& input){
  if(key == "rectilinear") return new PV_DiscreteRect(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}

AtomGroup* AtomGroup_Factory(std::string key, InputPack& input){
  if(key == "name") return new AtomGroup_name(input);
  if(key == "ndx") return new AtomGroup_ndx(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}

Calculation* Calculation_Factory(std::string key, InputPack& input){
  if(key == "isosurface") return new Calc_Isosurface(input);
  if(key == "nv") return new Calc_Nv(input);
  if(key == "angle") return new Calc_Angle(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}