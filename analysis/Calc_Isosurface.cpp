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
  Vec3<double> box_size;
  for(int i = 0; i < 3; i++){
    box_size[i] = box.boxvec[i][i];
  }
  if(frame_counter_ == 0){
    average_ = new VoxelGrid(npoints_, box_size, density_, sigma_, isovalue_);
    frame_ = new VoxelGrid(npoints_, box_size, density_, sigma_, isovalue_);
  }
  else{
    average_->setLength(box_size);
  }
  for(int i = 0; i < box.atoms.size(); i++ ){
    frame_->add_gaussian(box.atoms[i].x);
  }
  marchingCubes("golosio", *frame_, mesh_);
  return;
}
std::string Calc_Isosurface::printConsoleReport(){
  return;
}