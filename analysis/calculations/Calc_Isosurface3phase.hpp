#pragma once
#include "Calculation.hpp"
#include "../../cubes/MarchingCubesInterface.hpp"
class Calc_Isosurface3phase : public Calculation{
public:
  Calc_Isosurface3phase(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void output();
  virtual void update();
protected:
  virtual void printOutput();
  //number of actual frames computed
  int frame_counter_;
  bool initialized_ = 0;
  std::string method_;
  //total surface area of the isosurface
  std::vector<double> liquid_areas_, solid_areas_, sl_area, sv_area, lv_area, time_;
  //voxel grid input information
  Vec3<int> npoints_;
  Vec3<double> box_size_;
  double solid_area_, solid_sigma_, solid_density_, solid_isovalue_, distance_threshold_;
  double liquid_area_, liquid_sigma_, liquid_density_, liquid_isovalue_;
  //voxel grid output
  Mesh liquid_mesh_, solid_mesh_;
  VoxelGrid liquid_frame_, solid_frame_;
  VoxelGrid liquid_average_, solid_average_;
  AtomGroup* solid_atom_group_, *liquid_atom_group_; //these are the atoms for which calculations will be performed
  ProbeVolume* pv_; 
};
