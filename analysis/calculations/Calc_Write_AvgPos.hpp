//creates and writes to a gro file with the average atom positions and box positions for the specified atom group
//also reports some statistics about the system
#include "Calculation.hpp"
#include "../interface/interface.hpp" 
class Calc_Write_AvgPos : public Calculation{
public:
  Calc_Write_AvgPos(InputPack& input);
  virtual void update();
  virtual void calculate();
  virtual void output();
  virtual void finalOutput();

protected:
  void setRefPosArray();
  AtomGroup* atom_group_;
  int natoms_, frame_counter_;
  bool initialized_;
  std::vector<Vec3<double> > RefPos_, ComPos_; //contains initial atom position for computing pbcs
  std::vector<int> RefIdx_;
  std::string ofilename_;
  Vec3<double> box_size_, avg_box_size_;
};