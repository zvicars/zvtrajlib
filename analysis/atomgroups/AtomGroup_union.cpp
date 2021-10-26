#include "AtomGroup_union.hpp"
AtomGroup_union::AtomGroup_union(InputPack& input):AtomGroup{input}{
  std::string atom1name, atom2name;
  input.params().readString("atom_group1", KeyType::Required, atom1name);
  input.params().readString("atom_group2", KeyType::Required, atom2name);
  ag1_ = input.findAtomGroup(atom1name);
  ag2_ = input.findAtomGroup(atom2name);

  FANCY_ASSERT(ag1_ != 0, "Failed to find atomgroup 1 in AtomGroup_union");
  FANCY_ASSERT(ag2_ != 0, "Failed to find atomgroup 2 in AtomGroup_union");

  return;
}

void AtomGroup_union::update(){
  global_indices_.clear();
  auto gi1 = ag1_->getIndices();
  auto gi2 = ag2_->getIndices();
  global_indices_.insert(global_indices_.end(), gi1.begin(), gi1.end());
  global_indices_.insert(global_indices_.end(), gi2.begin(), gi2.end());
  return;
}