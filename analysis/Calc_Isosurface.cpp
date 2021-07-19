#include "Calc_Isosurface.hpp"
namespace CalculationRegistry {
static const Register<Calc_Isosurface>
  registerType("Isosurface");
}
Calc_Isosurface::Calc_Isosurface(InputPack& input):Calculation{input}
{
  std::vector<double> npoints;
  input.params().readVector("npoints", KeyType::Required, npoints);
  FANCY_ASSERT(npoints.size() == 3, "Invalid number of dimensions in npoints of isosurface calculation.");

  input.params().readNumber("sigma", KeyType::Required, sigma_);
  FANCY_ASSERT(sigma_ > 0, "Invalid sigma given for isosurface calculation.");

  input.params().readNumber("density", KeyType::Required, density_);
  FANCY_ASSERT(density_ > 0, "Invalid density given for isosurface calculation.");

  input.params().readNumber("isovalue", KeyType::Required, isovalue_);
  FANCY_ASSERT(isovalue_ > 0, "Invalid isovalue given for isosurface calculation.");  

  return;
}
Calc_Isosurface::~Calc_Isosurface(){
    delete frame_;
    delete average_;
    return;
}
void Calc_Isosurface::calculate(const Box& box){
  float sum = 0.0;
  for(int i = 0; i < box.atoms.size(); i++){
    sum += pv_->compute(box.atoms[i].x);
  }
  value_ = sum;
  return;
}
std::string Calc_Isosurface::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Nv" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}