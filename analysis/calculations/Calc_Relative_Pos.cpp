#include "Calc_Relative_Pos.hpp"
#include "../helper/pbc_functions.hpp"
#include "../../tools/stlmath.hpp"
Calc_Relative_Pos::Calc_Relative_Pos(InputPack& input) : Calculation_Histogram{input} {
  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  input.params().readString("grofile", KeyType::Required, refFile_);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  FANCY_ASSERT(!atom_group_->isDynamic(), "Calc_Relative_Pos requires the atom group to be static.");
  initialized_ = 0;
  frame_counter_ = 0;
  return;
}

void Calc_Relative_Pos::update(){
  if(hasUpdated()) return;
  Calculation::update();
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  if(!initialized_){
    natoms_ = atom_group_->getIndexCount();
    setRefPosArray(); //set initial atom positions to ensure correct pbc is used for edge atoms, also stores indices
    initialized_ = 1;
  }
  return;
}
void Calc_Relative_Pos::setRefPosArray(){
  Box gro_in;
  readGRO(refFile_, gro_in);
  refAtoms_ = gro_in.atoms;
  return;
}

void Calc_Relative_Pos::calculate(){
  if(!doCalculate()) return;
  std::vector<std::array<double,3> > dx_step(natoms_);
  auto& indices = atom_group_->getIndices();
  std::array<double,3> mean_step; mean_step.fill(0.0);
  std::array<double,3> var_step; var_step.fill(0.0);
  for(int i = 0; i < indices.size(); i++){
    int index = indices[i];
    auto newpos = getNearestPeriodicPosition(box->atoms[index].x, refAtoms_[index].x, box_size_);
    dx_step[i] = newpos - refAtoms_[index].x;
    mean_step = mean_step + dx_step[i];
    var_step = var_step + (dx_step[i] * dx_step[i]);
  }
  dx_vals_.insert(dx_vals_.begin(), dx_step.begin(), dx_step.end());
  mean_step = mean_step * (1.0/indices.size());
  var_step = var_step * (1.0/indices.size());
  means_.push_back(mean_step);
  vars_.push_back(var_step);
  times_.push_back(current_time_);
  frames_.push_back(current_frame_);
  frame_counter_++;
  return;
}

void Calc_Relative_Pos::output(){
  if(!doOutput()) return;
  return;
}

void Calc_Relative_Pos::finalOutput(){
  std::array<double,3> final_mean = {0.0,0.0,0.0}, final_var = {0.0,0.0,0.0};
  for(int i = 0; i < means_.size(); i++){
    final_mean = final_mean + means_[i];
    final_var = final_var + vars_[i];
  }
  final_mean = final_mean * (1.0/means_.size());
  final_var = final_var * (1.0/means_.size());

  std::vector<double> dx_vals;
  std::vector<double> dy_vals;
  std::vector<double> dz_vals;
  for(auto vals : dx_vals_){
    dx_vals.push_back(vals[0]);
    dy_vals.push_back(vals[1]);
    dz_vals.push_back(vals[2]);
  }
  std::vector<double> xx, xy, xz;
  std::vector<int> yx, yy, yz;
  makeHistogram(dx_vals, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, xx, yx);
  makeHistogram(dy_vals, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, xy, yy);
  makeHistogram(dz_vals, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, xz, yz);

  std::ofstream ofile(name_ + ".histogram");
  ofile << "Mean dx = " << final_mean[0] << "  " << final_mean[1] << "   " << final_mean[2] << "\n";
  ofile << "Var dx = " << final_var[0] << "  " << final_var[1] << "   " << final_var[2] << "\n";
  int nrows = std::max(xx.size(), std::max(xy.size(), xz.size()));
  for(int i = 0; i < nrows; i++){
    if(i < xx.size()) ofile << xx[i] << "   " << yx[i] << "   ";
    if(i < xy.size()) ofile << xy[i] << "   " << yy[i] << "   ";
    if(i < xz.size()) ofile << xz[i] << "   " << yz[i] << "   ";
    ofile << "\n";
  }
  ofile.close();
  return;
}