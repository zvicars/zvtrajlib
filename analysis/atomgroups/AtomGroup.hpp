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
  std::string getName() const {
    return name_;
  }
  std::string getType() const {
    return type_;
  }
  virtual void update(){
    return;
  }
protected:
  std::string name_;
  std::string type_;
  std::vector<int> global_indices_;
};