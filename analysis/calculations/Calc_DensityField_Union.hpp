#include "Calc_DensityField.hpp"
#include "../../interface/interface.hpp"
#include "../../tools/DifferentialGeometry.hpp"
#include "../../tools/stlmath.hpp"
#include "Calc_DensityField.hpp"
//does not inherit from densityfield because it doesn't need any of the book-keeping
//also needs to change some datatypes for subsequent analysis
class Calc_DensityField_Union : public Calculation{
  public:
  Calc_DensityField_Union(InputPack& input):Calculation{input}{
    input.params().readVector("names", KeyType::Required, field_names_);  
    for(auto name : field_names_){
      auto calc_pointer = input.findCalculation(name);
      FANCY_ASSERT(calc_pointer != 0, "Failed to find specified calculation.");
      //static cast should fail if it isn't the appropriate derived class
      Calc_DensityField* df_ptr = static_cast<Calc_DensityField*>(calc_pointer);
      df_vec_.push_back(df_ptr);
    }
    std::vector<std::size_t> npoints;
    Vec3<std::size_t> npoints_arr;
    input.params().readVector("npoints", KeyType::Required, npoints);
    FANCY_ASSERT(npoints.size() == 3, "Inappropriate npoints specified.");
    for(std::size_t i = 0; i < 3; i++){
      npoints_arr[i] = npoints[i];
    }
    data_matrix_.initialize(npoints_arr);
    return;
  }
  virtual void finalOutput();
  protected:
  Matrix<double,3> data_matrix_;
  std::vector<Calc_DensityField*> df_vec_;
  std::vector<std::string> field_names_;
};