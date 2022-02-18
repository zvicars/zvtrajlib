#include "AtomGroup_intersection.hpp"
AtomGroup_intersection::AtomGroup_intersection(InputPack& input):AtomGroup{input}{
  std::string atom1name, atom2name;
  input.params().readString("atom_group1", KeyType::Required, atom1name);
  input.params().readString("atom_group2", KeyType::Required, atom2name);
  ag1_ = input.findAtomGroup(atom1name);
  ag2_ = input.findAtomGroup(atom2name);

  FANCY_ASSERT(ag1_ != 0, "Failed to find atomgroup 1 in AtomGroup_intersection");
  FANCY_ASSERT(ag2_ != 0, "Failed to find atomgroup 2 in AtomGroup_intersection");

  return;
}

void AtomGroup_intersection::update(){
  if(hasUpdated()) return;
  AtomGroup::update();
  if(!ag1_->hasUpdated()) ag1_->update();
  if(!ag2_->hasUpdated()) ag2_->update();  
  index_set_.clear();
  global_indices_.clear();
  auto gi1 = ag1_->getIndices();
  auto gi2 = ag2_->getIndices();
  for(int i = 0; i < gi1.size(); i++){
    index_set_.insert(gi1[i]);
  }
  for(int i = 0; i < gi2.size(); i++){
    //see if the key already exists
    auto result = index_set_.find(gi2[i]);
    //if it does, we know it's in both atomgroups
    if(result != index_set_.end()) global_indices_.push_back(gi2[i]);
  }
  return;
}