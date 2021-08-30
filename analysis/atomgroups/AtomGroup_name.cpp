#include "AtomGroup_name.hpp"
AtomGroup_name::AtomGroup_name(InputPack& input):AtomGroup{input}{
  input.params().readString("atom_name", KeyType::Required, atomname_);
  auto box = input.getBox();
  FANCY_ASSERT(box->hasNamedAtoms == 1, "Box atoms are not named, name-type atomgroups cannot be used unless you provide a gro/top file.");
  for(std::size_t i = 0; i < box->atoms.size(); i++){
    if(box->atoms[i].name == atomname_){
      global_indices_.push_back(i-1);
    }
  }
  if(global_indices_.size() == 0){
    std::cout << "Warning... this atomgroup does not have any atoms in it!" << std::endl;
  }
  return;
}