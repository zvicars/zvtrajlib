#pragma once
#include "Calculation.hpp"
#include "../../cubes/MarchingCubesInterface.hpp"
class Calc_IsosurfaceMultiphase : public Calculation{
public:
  Calc_IsosurfaceMultiphase(InputPack& input);
  ~Calc_IsosurfaceMultiphase(){
    if(!pv_specified_) delete pv_; 
  }
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void output();
  virtual void update();
protected:
  virtual void printOutput();
  //number of actual frames computed
  int frame_counter_, num_groups_;
  bool initialized_ = 0, pv_specified_=0;
  std::string method_;
  //total surface area of the isosurface
  std::vector<std::vector<double> >  area_vectors_;
  //voxel grid input information
  Vec3<int> npoints_;
  Vec3<double> box_size_;
  std::vector<double> areas_, sigmas_, densities_, isovalues_, distance_rmax_, distance_sigmas_, num_neighbor_thresholds_, num_neighbor_sigmas_;
  //voxel grid output
  std::vector<std::string> filepaths_;
  std::vector<Mesh> meshes_;
  std::vector<VoxelGrid> frames_;
  std::vector<AtomGroup*> atomgroups_; //these are the atoms for which calculations will be performed
  std::vector<std::string> agnames_, hex_colors_;
  bool outputPLY_ = 0;
  //ts output
  std::vector<double> times_;
  std::vector<std::vector<double> > avecs_;


  ProbeVolume* pv_; 
};
