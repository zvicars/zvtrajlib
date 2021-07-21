#include "Calc_Isosurface.hpp"
namespace CalculationRegistry {
static const Register<Calc_Isosurface>
  registerType("Isosurface");
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
  input.params().readNumber("sigma", KeyType::Required, sigma_);
  FANCY_ASSERT(sigma_ > 0, "Invalid sigma given for isosurface calculation.");

  input.params().readNumber("density", KeyType::Required, density_);
  FANCY_ASSERT(density_ > 0, "Invalid density given for isosurface calculation.");

  input.params().readNumber("isovalue", KeyType::Required, isovalue_);
  FANCY_ASSERT(isovalue_ > 0, "Invalid isovalue given for isosurface calculation.");  

  return;
}
void Calc_Isosurface::calculate(const Box& box){
  current_time_ = box.time;
  current_frame_ = box.frame_counter;
  if(output_freq_ <= 0) return;
  if(current_time_ < equilibration_) return;
  if(current_frame_%output_freq_ != 0) return;
  Vec3<double> box_size;
  for(int i = 0; i < 3; i++){
    box_size[i] = box.boxvec[i][i];
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
  for(int i = 0; i < box.atoms.size(); i++ ){
    frame_.add_gaussian(box.atoms[i].x);
  }
  marchingCubes("golosio", frame_, mesh_);
  //average_->sumInPlace(frame_);
  printOutput();
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