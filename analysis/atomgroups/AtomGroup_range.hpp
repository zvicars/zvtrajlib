#pragma once
#include "AtomGroup.hpp"

class AtomGroup_range : public AtomGroup{
public:
  AtomGroup_range(InputPack& input);
private:
  std::vector<double> params_;
};