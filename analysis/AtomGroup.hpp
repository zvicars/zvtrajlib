#pragma once
#include <vector>
#include <string>
#include "../tools/GenericFactory.hpp"
#include "InputPack.hpp"

class AtomGroup{
public:
  using KeyType = ParameterPack::KeyType;
  AtomGroup(InputPack&);
  const std::vector<int>& getIndices(){
    return global_indices_;
  }
  std::string getName(){
    return name_;
  }
protected:
  std::string name_;
  std::string type_;
  std::vector<int> global_indices_;
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
