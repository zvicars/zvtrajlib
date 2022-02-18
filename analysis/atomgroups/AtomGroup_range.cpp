#include "AtomGroup_range.hpp"
AtomGroup_range::AtomGroup_range(InputPack& input):AtomGroup{input}{
  input.params().readVector("params", KeyType::Required, params_);
  auto box = input.getBox();
  for(std::size_t i = params_[0]; i < params_[1]; i+= params_[2]){
    int idx = i-1;
    if(idx >= box->atoms.size()) continue;
    global_indices_.push_back(idx);
  }
  if(global_indices_.size() == 0){
    std::cout << "Warning... this atomgroup does not have any atoms in it!" << std::endl;
  }
  return;
}