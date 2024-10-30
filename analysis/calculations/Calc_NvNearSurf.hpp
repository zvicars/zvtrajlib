#pragma once
#include "Calc_Histogram.hpp"
#include "../../tools/cellgrid.hpp"
#include <cmath>
//uses ice-id style arguments to count the total number of ice-like waters excluding those that absolutely cannot be ice
class Calc_NvNearSurf : public Calculation_Histogram{
public:
  Calc_NvNearSurf(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void output(){
    Calculation_Histogram::output();
    if(doOutput()){
      if(bOutput) printXYZ();
    }
  }
  virtual void finalOutput();
  virtual void update(){
    Calculation::update();
    for(int i = 0; i < 3; i++){
      box_size_[i] = box->boxvec[i][i];
    }
    return;
  }
protected:
  void updateNeighborLists(){
    const auto& surf_atom_indices = surf_group_->getIndices();
    const auto& atom_positions = box->atoms;
    c_surf.reset(surface_radius_, box_size_);
    for(auto index : surf_atom_indices){
      c_surf.addIndexToGrid(index, atom_positions[index].x);
    }

    const auto& ice_atom_indices = ice_group_->getIndices();
    c_ice.reset(ice_radius_, box_size_);
    for(auto index : ice_atom_indices){
      c_ice.addIndexToGrid(index, atom_positions[index].x);
    }

    return;
  }
  double neighborEval(const Vec3<double>& pos, CellGrid& cg, double radius, double smear, double cutoff) const; 
  void printXYZ();
  CellGrid c_surf, c_ice;
  Vec3<double> box_size_;
  Vec<double> time_vec_;
  Vec<int> step_vec_;
  Vec<double> count_vec_;
  double value_, mean_, var_;
  std::string pv_name_;
  ProbeVolume* pv_;
  AtomGroup *atom_group_, *surf_group_, *ice_group_; //these are the atoms for which calculations will be performed
  int dump_frame_ = -1;

  double ice_radius_, surface_radius_;
  double surface_thresh_, bulk_ice_thresh_, bridging_ice_thresh_;
  double ice_smear_, surface_smear_;

  //atom-wise data vectors
  Vec<double> ice_neighbors_, surf_neighbors_, pv_eval_;

  //printing ouput files
  std::string trajfile_;
  std::ofstream trajout_;
  bool bOutput=0;

};