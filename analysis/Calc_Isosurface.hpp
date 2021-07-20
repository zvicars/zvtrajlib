#pragma once
#include "Calculation.hpp"
#include "../cubes/MarchingCubesInterface.hpp"
class Calc_Isosurface : public Calculation{
public:
  Calc_Isosurface(InputPack& input);
  ~Calc_Isosurface(){
      delete frame_;
      delete average_;
  }
  virtual void calculate(const Box& box);
  virtual std::string printConsoleReport();
  virtual void printOutput();
  virtual void finalOutput();
private:
  //number of actual frames computed
  int frame_counter_;
  //total surface area of the isosurface
  std::vector<double> areas_;
  //voxel grid input information
  Vec3<int> npoints_;
  double area_, sigma_, density_, isovalue_;
  //voxel grid output
  Mesh mesh_;
  VoxelGrid* frame_;
  VoxelGrid* average_;
};
