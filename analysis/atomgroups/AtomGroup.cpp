#include "AtomGroup.hpp"
AtomGroup::AtomGroup(InputPack& input){
  input.params().readString("name", KeyType::Required, name_);
  input.params().readString("type", KeyType::Required, type_);
  input.addAtomGroup(name_, this);  
  return;
}