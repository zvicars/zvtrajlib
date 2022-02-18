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
  int getIndexCount() const {
    return global_indices_.size();
  }
  int getIndex(int i) const {
    return global_indices_[i];
  }
  std::string getName() const {
    return name_;
  }
  std::string getType() const {
    return type_;
  }
  virtual void update(){
    update_flag_ = 1;
    return;
  }
  void finish(){
    update_flag_ = 0;
  }
  virtual bool isDynamic() const {
    return 0;
  }
  bool hasUpdated() const {
    return update_flag_;
  }
protected:
  std::string name_;
  std::string type_;
  std::vector<int> global_indices_;
  bool update_flag_;
};