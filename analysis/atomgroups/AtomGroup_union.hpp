#pragma once
#include "AtomGroup.hpp"
#include <unordered_set>

class AtomGroup_union : public AtomGroup{
public:
  AtomGroup_union(InputPack& input);
  virtual void update();
private:
  AtomGroup* ag1_;
  AtomGroup* ag2_;
};