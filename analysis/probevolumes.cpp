#include "probevolumes.hpp"
//place in registry
ProbeVolume* registerPV_DiscreteRect(std::string input){
  return new PV_DiscreteRect(std::string input);
}
reg::probe_volumes.registerParser( "PV_DiscreteRect", &registerPV_DiscreteRect );

class ProbeVolume{ //probe volumes allow you to specify a particle position and receive a floating point number corresponding to whether the particle is in the volume or not
public:
  ProbeVolume(std::string input){
    viaKey<std::string>("name", input, name_);
    return;
  }
  virtual double compute(Vec3<double> position) = 0;
private:
  std::string name_;
};

class PV_DiscreteRect : public ProbeVolume{
public: 
  PV_DiscreteRect(std:string input):ProbeVolume{input}
  {
    Vec<double> corners;
    viaKeyVec<double>("corners", input, corners);
    if(corners.size() != 6){
      std::cout << "Invalid number of arguments given for PV_DiscreteRect parameter \"corners\"." << std::endl;
      throw 0;
    }
    for(int i = 0; i < 3; i++){
      lower_corner_[i] = corners[i];
      upper_corner_[i] = corners[i+3];
    }
  }
  double calculate(Vec3<double> position){
    for(int i = 0; i < 3; i++){
      if( position[i] < lower_corner_[i] || position[i] > upper_corner_[i] ) return 0.0; 
    }
    return 1.0;
  }
private:
  Vec3<double> lower_corner_;
  Vec3<double> upper_corner_;
};
