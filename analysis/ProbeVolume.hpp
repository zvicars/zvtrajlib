#pragma once
#include "../tools/Assert.hpp"
#include "../tools/GenericFactory.hpp"
#include "../interface/datatypes.hpp"
#include "InputPack.hpp"

class ProbeVolume{ //probe volumes allow you to specify a particle position and receive a floating point number corresponding to whether the particle is in the volume or not
public:
  using KeyType = ParameterPack::KeyType;
  ProbeVolume(InputPack& input);
  virtual ~ProbeVolume(){
    return;
  }
  virtual double compute(Vec3<double> position) = 0;
private:
  std::string name_;
};

namespace ProbeVolumeRegistry {

// Manages the creation of Steinhardt objects
using Factory = GenericFactory<
  ProbeVolume,   // base class
  std::string,  // key type
  InputPack&  // input types
>;

// Object that registers the mapping with the factory
template<typename S>
using Register = RegisterInFactory<
  ProbeVolume,   // base class
  S,            // derived class
  std::string,  // key type
  InputPack&
>;
}
