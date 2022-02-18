#pragma once
#include "AtomGroup.hpp"
#include <unordered_set>

class AtomGroup_union : public AtomGroup{
public:
  AtomGroup_union(InputPack& input);
  virtual void update();
  virtual bool isDynamic(){
    for(auto ag : ags_){
      if(ag->isDynamic()) return 1;
    }
    return 0;
  }
private:
  std::unordered_set<int> index_set_;
  std::vector<AtomGroup*> ags_;
};