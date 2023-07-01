#pragma once
#include "Calc_Histogram.hpp"
#include "op_components/qlmi.hpp"
#include "op_components/htildei.hpp"
#include "op_components/switch.hpp"
#include "../../tools/cellgrid.hpp"
class Calc_Mv : public Calculation_Histogram{
public:
  Calc_Mv(InputPack& input);
  ~Calc_Mv(){
    delete probe_volume_;
    delete qlmi_;
    delete qbar_cutoff_;
    delete h4nn_switch_;
    delete shell_;
  }
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void getAtomwiseInformation(Vec<Vec3<double> >& positions, Vec<double>& mvals){
    positions_ = positions;
    mvals_ = mvals;
    return;
  }
  virtual void update(){
    Calculation::update();
    return;
  }
  void setAtomwise(){
    atomwiseNeeded = 1;
    return;
  }
protected:
  void buildCellGrid(const std::vector<Vec3<double> >& positions){
    c1_.reset(qlmi_rmax_ + qlmi_rcut_, box_size_);
    for(int i = 0; i < positions.size(); i++){
      c1_.addIndexToGrid(i, positions[i]);
    }
  }
  void getNeighbors(int idx, Vec3<double> position, Vec<int>& neighbors){
    neighbors.clear();
    neighbors = c1_.getNearbyIndices(position);
    for(int i = 0; i < neighbors.size(); i++){
      if(neighbors[i] == idx) neighbors.erase(neighbors.begin()+i);
    }
    return;
  }
  Vec3<double> box_size_;
  Vec<double> time_vec_;
  Vec<int> step_vec_;
  Vec<double> count_vec_;
  double value_, mean_, var_;
  std::string pv_name_;
  AtomGroup* atom_group_; //these are the atoms for which calculations will be performed
  int dump_frame_ = -1;
  simpleTwoSidedSwitch* h4nn_switch_;
  simpleSwitch* qbar_cutoff_;
  cuboidalINDUSVolume* probe_volume_, * shell_;
  qlmi* qlmi_;
  double qlmi_rmax_,qlmi_rcut_, qlmi_order_; //for testing
  bool atomwiseNeeded = 0;
  Vec<Vec3<double> > positions_;
  Vec<double> mvals_;
  CellGrid c1_;
};