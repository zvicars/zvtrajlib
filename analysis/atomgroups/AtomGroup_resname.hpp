#pragma once
#include "AtomGroup.hpp"

class AtomGroup_resname : public AtomGroup{
public:
  AtomGroup_resname(InputPack& input);
private:
  std::string resname_;
};