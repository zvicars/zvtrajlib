#pragma once
#include "Calculation.hpp"
#include "../helper/make_histogram.hpp"
#include <cmath>
class Calc_Lattice_Motion : public Calculation{
public:
  Calc_Lattice_Motion(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void update();
  double getNearestPeriodicImage(Vec3<double>& pos, const Vec3<double>& ref_pos);
private:
  //logic
  bool isInitialized; //has the first frame been loaded?
  //particle identification
  std::string atom_group_name_;
  AtomGroup* atom_group_;
  int ref_atom_index_;
  Vec3<double> ref_atom_position_; //set by update()
  
  int num_nn_;
  Vec<int> nn_indices_; //set by update()
  Vec< Vec3<double> > nn_positions_; //set by update()
  Vec<double> nn_distances_;
  Vec3<double> nn_original_com_;
  //statistical information
  Vec<Vec3<double> > offsets_;
  Vec<int> step_vec_;
  Vec<double> time_vec_;
  //misc stuff for calculations
  Vec3<double> box_size_;
};