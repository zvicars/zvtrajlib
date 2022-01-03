#pragma once
#include "AtomGroup.hpp"
#include <unordered_set>

class AtomGroup_not : public AtomGroup{
public:
  AtomGroup_not(InputPack& input);
  virtual void update();
  virtual bool isDynamic(){
    return ag1_->isDynamic() || ag2_->isDynamic();
  }
private:
  std::unordered_set<int> index_set_;
  AtomGroup* ag1_;
  AtomGroup* ag2_;
};