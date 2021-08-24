#include "Calc_SWIPES_CircleFit.hpp"
/*
//fits a circle to a 2d grid density
class Calc_CircleFit : public Calculation{
public:
  Calc_CircleFit(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:
  Vec<Vec<double> > angle_distribution_;
  Vec<double> angle_vec_;
  Vec<int> step_vec_;
  Vec<double> time_vec_;
  double mean_, var_, value_;
  //this one depeends on 3 atomgroups of equal size
  Vec3<AtomGroup*> atom_groups_;
};*/


Calc_SWIPES_CircleFit::Calc_SWIPES_CircleFit(InputPack& input) : Calculation{input} {
  //boudary within which to consider atoms
  std::string pv_name_;
  input.params().readString("probe_volume", KeyType::Required, pv_name_);
  auto pv_pointer = input.findProbeVolume(pv_name_);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer;

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");

  Vec<double> xmin, xmax;
  input.params().readVector("xmin", KeyType::Required, xmin);
  input.params().readVector("xmax", KeyType::Required, xmax);
  FANCY_ASSERT(xmin.size() == 2, "Invalid xmin size in CircleFit object.");
  FANCY_ASSERT(xmax.size() == 2, "Invalid xmax size in CircleFit object.");
  for(int i = 0; i < 2; i++){
    xmin_[i] = xmin[i];
    xmax_[i] = xmax[i];
  }

  std::vector<double> params(3, 1.0);
  input.params().readVector("params", KeyType::Optional, params);
  FANCY_ASSERT(params.size() == 3, "Invalid params size in CircleFit object.");
  for(int i = 0; i < 3; i++){
    params_[i] = params[i];
  }

  std::vector<bool> fix_flags(3, 0);
  input.params().readVector("fix_params", KeyType::Optional, fix_flags);
  FANCY_ASSERT(fix_flags.size() == 3, "Invalid fix_params size in CircleFit object.");
  for(int i = 0; i < 3; i++){
    fix_params_[i] = fix_flags[i];
  } 
  input.params().readNumber("density", KeyType::Optional, density_);
  return;

}

void Calc_2D_Density::finalOutput(){
  auto grid_data = calc_;
  return;
}