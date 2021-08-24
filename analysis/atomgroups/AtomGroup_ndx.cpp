#include "AtomGroup_ndx.hpp"
AtomGroup_ndx::AtomGroup_ndx(InputPack& input):AtomGroup{input}{
  input.params().readString("label", KeyType::Required, label_);
  auto box = input.getBox();
  auto it = box->idxinfo.indexes.find(label_);
  FANCY_ASSERT(it != box->idxinfo.indexes.end(), "Failed to find a specified index in the indexfile.")
  global_indices_ = it->second;
  for(int i = 0; i < global_indices_.size(); i++){
    global_indices_[i] -= 1;
  }
  return;
}