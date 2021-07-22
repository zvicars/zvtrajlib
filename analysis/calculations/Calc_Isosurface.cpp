#include "Calc_Isosurface.hpp"
namespace CalculationRegistry {
static const Register<Calc_Isosurface>
  registerType("isosurface");
}

Calc_Isosurface::Calc_Isosurface(InputPack& input):Calculation{input}
{
  frame_counter_ = 0;
  std::vector<double> npoints;
  input.params().readVector("npoints", KeyType::Required, npoints);
  FANCY_ASSERT(npoints.size() == 3, "Invalid number of dimensions in npoints of isosurface calculation.");
  for(int i = 0; i < 3; i++){
    npoints_[i] = npoints[i];
  }

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  
  input.params().readNumber("sigma", KeyType::Required, sigma_);
  FANCY_ASSERT(sigma_ > 0, "Invalid sigma given for isosurface calculation.");
  input.params().readNumber("density", KeyType::Required, density_);
  FANCY_ASSERT(density_ > 0, "Invalid density given for isosurface calculation.");
  input.params().readNumber("isovalue", KeyType::Required, isovalue_);
  FANCY_ASSERT(isovalue_ > 0, "Invalid isovalue given for isosurface calculation.");
  method_ = "golosio";
  input.params().readString("method", KeyType::Optional, method_);
  FANCY_ASSERT(method_ == "golosio", "Invalid method chosen for instantaneous interface calculation, valid options are \'golosio\' and \'rchandra\'.");
  return;
}
void Calc_Isosurface::calculate(){
  current_time_ = box->time;
  current_frame_ = box->frame_counter;
  if(!doCalculate()) return;
  Vec3<double> box_size;
  for(int i = 0; i < 3; i++){
    box_size[i] = box->boxvec[i][i];
  }
  if(!initialized_){
    average_.initialize(npoints_, box_size, density_, sigma_, isovalue_);
    frame_.initialize(npoints_, box_size, density_, sigma_, isovalue_);
    initialized_ = 1;
  }
  else{
    frame_.setLength(box_size);
    frame_.clear();
  }
  for(int i = 0; i < atom_group_->getIndices().size(); i++ ){
    int idx = atom_group_->getIndices()[i];
    frame_.add_gaussian(box->atoms[idx].x);
  }
  marchingCubes(method_, frame_, mesh_);
  //average_->sumInPlace(frame_);
  if(doOutput()) printOutput();
  frame_counter_++;
  return;
}

std::string Calc_Isosurface::printConsoleReport(){
  return "";
}

void Calc_Isosurface::printOutput(){
  std::string filepath = base_ + "_frame" + std::to_string(frame_counter_) + ".stl";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface step calculation.");
  std::string output;
  printSTL(mesh_, output);
  ofile << output;
  ofile.close();
};
void Calc_Isosurface::finalOutput(){
  return;
  if(output_freq_ <= 0) return;
  average_.scalarMult(1.0/(double)frame_counter_);
  std::string filepath = base_ + "_average.stl";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface average.");
  std::string output;
  printSTL(mesh_, output);
  ofile << output;
  ofile.close(); 
  return;
};