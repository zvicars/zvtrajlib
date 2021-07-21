#pragma once
#include <vector>
#include <string>
#include "../tools/GenericFactory.hpp"
#include "InputPack.hpp"

class AtomGroup{
public:
  AtomGroup(InputPack&);
  const std::vector<int>& getIndices();
protected:
  std::string name;
  std::string type;
  std::vector<int> global_indices;
};
namespace AtomGroupRegistry {
// Manages the creation of Steinhardt objects
using Factory = GenericFactory<
  AtomGroup,   // base class
  std::string,  // key type
  InputPack&  // input types
>;

// Object that registers the mapping with the factory
template<typename S>
using Register = RegisterInFactory<
  AtomGroup,   // base class
  S,            // derived class
  std::string,  // key type
  InputPack&
>;
}
