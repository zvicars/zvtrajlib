//creates and writes to an xtc file with all atoms at -1, -1, -1 except for ones in the atomgroup
//a slightly easier method of rendering solid-like atoms than using Sean's scripts
#pragma once
#include "Calculation.hpp"
#include "../../interface/interface.hpp"
class Calc_Write_Xtc : public Calculation{
public:
  Calc_Write_Xtc(InputPack& input);
  virtual void update();
  virtual void output();
  virtual void finalOutput();
protected:
  AtomGroup* atom_group_;
  xdr::XDRFILE* output_handle_;
  int xdr_natoms_, xdr_step_;
  float xdr_time_;
  xdr::matrix xdr_box_;
  xdr::rvec* xdr_x_;
  float xdr_prec_;
  bool initialized_;
  std::string filename_;
};