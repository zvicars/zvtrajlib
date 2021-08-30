#pragma once
#include "AtomGroup.hpp"
#include <unordered_set>

class AtomGroup_intersection : public AtomGroup{
public:
  AtomGroup_intersection(InputPack& input);
  virtual void update();
private:
  std::unordered_set<int> index_set_;
  AtomGroup* ag1_;
  AtomGroup* ag2_;
};