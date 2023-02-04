#pragma once
#include "../../tools/Assert.hpp"
#include "../../interface/datatypes.hpp"
#include "../InputPack.hpp"
class ProbeVolume{ //probe volumes allow you to specify a particle position and receive a floating point number corresponding to whether the particle is in the volume or not
public:
  using KeyType = ParameterPack::KeyType;
  ProbeVolume(InputPack& input);
  ProbeVolume(std::string name);
  ProbeVolume();
  virtual ~ProbeVolume(){
    return; 
  }
  virtual double compute(Vec3<double> position) = 0;
  virtual void update(){
    update_flag_ = 1;
    return;
  }
  void finish(){
    update_flag_ = 0;
  }
  bool hasUpdated(){
    return update_flag_;
  }
protected:
  const Box* box;
  std::string name_;
  bool update_flag_;
};