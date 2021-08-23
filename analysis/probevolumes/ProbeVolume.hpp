#pragma once
#include "../../tools/Assert.hpp"
#include "../../interface/datatypes.hpp"
#include "../InputPack.hpp"

class ProbeVolume{ //probe volumes allow you to specify a particle position and receive a floating point number corresponding to whether the particle is in the volume or not
public:
  using KeyType = ParameterPack::KeyType;
  ProbeVolume(InputPack& input);
  virtual ~ProbeVolume(){
    return; 
  }
  virtual double compute(Vec3<double> position) = 0;
  virtual void update(){
    return;
  }
protected:
  const Box* box;
private:
  std::string name_;
};