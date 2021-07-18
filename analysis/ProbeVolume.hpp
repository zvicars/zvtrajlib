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



// Register ProbeVolumes using a GenericFactory
// - Use these typedefs for brevity
//Generic factory implementation shamelessly borrowed from Sean's OP code
namespace ProbeVolumeRegistry {
template<typename P>
using Register = RegisterInFactory<
  ProbeVolume,            // base class
  P,                      // derived class
  std::string,            // generating key
  const InputPack&   // input type(s)
>;

using Factory = GenericFactory< ProbeVolume, std::string, const InputPack& >;
} // end namespace ProbeVolumeRegistry