//AtomGroup Dynamic: Accepts a timeseries of indices where each index is supposed to correspond to a particular frame
//this is useful for something like modeling ice or another dynamic quantity
//format is time (ps) idx_1 idx_2 ... idx_n so the script will do its best to match the current simulation time with the nearest index
#pragma once
#include "AtomGroup.hpp"
class AtomGroup_dynamic : public AtomGroup{
public:
  AtomGroup_dynamic(InputPack&);
  virtual void update();
protected:
  std::string filename_;
  std::vector<std::vector<int> > frames_;
  std::vector<double> times_;
  double dt_, tmin_, tmax_;
  const Box* box_; //this one needs the box pointer
  bool allowOB_;
  int index_;
};