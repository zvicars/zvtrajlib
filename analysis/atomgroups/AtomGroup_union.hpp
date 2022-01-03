#pragma once
#include "AtomGroup.hpp"
#include <unordered_set>

class AtomGroup_union : public AtomGroup{
public:
  AtomGroup_union(InputPack& input);
  virtual void update();
  virtual bool isDynamic(){
    return ag1_->isDynamic() || ag2_->isDynamic();
  }
private:
  AtomGroup* ag1_;
  AtomGroup* ag2_;
};