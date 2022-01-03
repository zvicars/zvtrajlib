#include "AtomGroup_system.hpp"
AtomGroup_system::AtomGroup_system(InputPack& input):AtomGroup{input}{
  auto box = input.getBox();
  for(std::size_t i = 0; i < box->atoms.size(); i++){
    global_indices_.push_back(i);
  }
  if(global_indices_.size() == 0){
    std::cout << "Warning... this atomgroup does not have any atoms in it!" << std::endl;
  }
  return;
}