#pragma once
#include "../tools/Assert.hpp"
#include "../tools/GenericFactory.hpp"
#include "InputPack.hpp"
#include "../interface/datatypes.hpp"

class ProbeVolume{ //probe volumes allow you to specify a particle position and receive a floating point number corresponding to whether the particle is in the volume or not
public:
  using KeyType = ParameterPack::KeyType;
  ProbeVolume(InputPack& input);
  virtual double compute(Vec3<double> position) = 0;
private:
  std::string name_;
};

namespace ProbeVolumeRegistry {

using Key  = std::string;
using Base = ProbeVolume;
using Factory = GenericFactory<Base, Key, InputPack&>;


template<typename D>
using Register = RegisterInFactory<Base, D, Key, InputPack&>;
}