//Atomgroup corresponding to loaded ndx entries, requires an ndx file to have been loaded
#pragma once
#include "AtomGroup.hpp"

class AtomGroup_ndx : public AtomGroup{
public:
  AtomGroup_ndx(InputPack& input);
private:
  std::string label_;
};