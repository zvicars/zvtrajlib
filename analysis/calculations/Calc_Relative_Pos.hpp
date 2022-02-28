//creates and writes to a gro file with the average atom positions and box positions for the specified atom group
//also reports some statistics about the system
#pragma once
#include "Calc_Histogram.hpp"
#include "../interface/interface.hpp" 
class Calc_Relative_Pos: public Calculation_Histogram{
public:
  Calc_Relative_Pos(InputPack& input);
  virtual void update();
  virtual void calculate();
  virtual void output();
  virtual void finalOutput();

protected:
  void setRefPosArray();
  AtomGroup* atom_group_;
  int natoms_, frame_counter_;
  bool initialized_;
  std::string refFile_;
  std::vector<Atom> refAtoms_;
  Vec3<double> box_size_;
  std::vector<double> times_;
  std::vector<double> frames_;
  std::vector<std::array<double, 3> > means_, vars_, dx_vals_;
};