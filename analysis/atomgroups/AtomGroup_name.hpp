#pragma once
#include "../AtomGroup.hpp"

class AtomGroup_name : public AtomGroup{
public:
  AtomGroup_name(InputPack& input);
private:
  std::string atomname_;
};