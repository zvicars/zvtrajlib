//uses a fit on the 1D density profile to set an additional boundary to a calcNv calculation
#pragma once
#include "Calc_Nv.hpp"
#include "Calc_1D_Density.hpp"
class Calc_Nv_wFit : public Calc_Nv{
public:
  Calc_Nv_wFit(InputPack& InputPack);
  virtual void calculate();
protected:
  Calc_1D_Density* calc_;
  bool dir_; //0 is less than, 1 is greater than
};