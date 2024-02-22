#include "PV_Trivial.hpp"
PV_Trivial::PV_Trivial(InputPack& input):ProbeVolume{input}
{
  return;
}

PV_Trivial::PV_Trivial(){
  return;
}

double PV_Trivial::compute(Vec3<double> position){
  return 1.0;
}