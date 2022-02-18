#include "AtomGroup_union.hpp"
AtomGroup_union::AtomGroup_union(InputPack& input):AtomGroup{input}{
  std::vector<std::string> agnames;
  input.params().readVector("atom_groups", KeyType::Required, agnames);
  for(auto name : agnames){
    auto ag = input.findAtomGroup(name);
    FANCY_ASSERT(ag != 0, "Failed to find atomgroup " + name + " in AtomGroup_union");
    ags_.push_back(ag);
  }
  return;
}

void AtomGroup_union::update(){
  if(hasUpdated()) return;
  AtomGroup::update();
  for(auto ag : ags_){
    if(!ag->hasUpdated()) ag->update();
  }
  index_set_.clear();
  global_indices_.clear();
  for(auto ag : ags_){
    auto gi = ag->getIndices();
    for(auto& index : gi){
      index_set_.insert(index);
    }
  }
  for(auto& index : index_set_){
    global_indices_.push_back(index);
  }
  return;
}