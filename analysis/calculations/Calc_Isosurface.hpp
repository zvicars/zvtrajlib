#pragma once
#include "Calculation.hpp"
#include "../../cubes/MarchingCubesInterface.hpp"
class Calc_Isosurface : public Calculation{
public:
  Calc_Isosurface(InputPack& input);
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
  std::vector<double> areas_;
  //voxel grid input information
  Vec3<int> npoints_;
  double area_, sigma_, density_, isovalue_;
  //voxel grid output
  Mesh mesh_;
  VoxelGrid frame_;
  VoxelGrid average_;
  AtomGroup* atom_group_; //these are the atoms for which calculations will be performed

};
