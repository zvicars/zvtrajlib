#pragma once
#include <vector>
#include <string>
#include "../InputPack.hpp"

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