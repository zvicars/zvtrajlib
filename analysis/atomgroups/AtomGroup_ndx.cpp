#include "AtomGroup_ndx.hpp"
namespace AtomGroupRegistry {
static const Register<AtomGroup_ndx>
  registerType("ndx");
}
AtomGroup_ndx::AtomGroup_ndx(InputPack& input):AtomGroup{input}{
  input.params().readString("label", KeyType::Required, label_);
  box = input.getBox();
  auto it = box->idxinfo.indexes.find(label_);
  FANCY_ASSERT(it != box->idxinfo.indexes.end(), "Failed to find a specified index in the indexfile.")
  global_indices_ = it->second;
  return;
}