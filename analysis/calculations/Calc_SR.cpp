#include "Calc_SR.hpp"
#include "../tools/stlmath.hpp"
Calc_SR::Calc_SR(InputPack& input):Calculation_Histogram{input}{
  input.params().readNumber("cutoff", KeyType::Required, cutoff_);
  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group."); 
  U_step_.resize(box->atoms.size()); 
  return;
}
void Calc_SR::update(){
  Calculation_Histogram::update();
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  return;
}
void Calc_SR::calculate(){
  if(!doCalculate()) return;
  calculate_flag_ = 1;
  auto indices = atom_group_->getIndices();
  //update cell grid
  CellGrid c1;
  c1.reset(cutoff_, box_size_);
  for(auto idx : indices){
    c1.addIndexToGrid(idx, box->atoms[idx].x);
  }
  std::vector<bool> hasComputed(box->atoms.size(), 0);
  for(auto& val : U_step_){
    val = 0.0;
  }
  //run through each pair of atoms and calculate potential between them
  for(auto idx : indices){
    auto indices2 = c1.getNearbyIndices(idx, box->atoms[idx].x);
    for(auto idx2 : indices2){
      if(hasComputed[idx2]) continue;
      double eval_step = compute_potential(idx, idx2);
      U_step_[idx] += eval_step;
      U_step_[idx2] += eval_step;
    }
    hasComputed[idx] = 1;
  }
  value_ = 0.0;
  for(auto val :  U_step_){
    value_ += val;
  }
  value_ *= 0.5;
  value_vec_.push_back(value_);
  time_vec_.push_back(box->time);
  step_vec_.push_back(box->frame);
  return;
}
void Calc_SR::finalOutput(){
  if(output_freq_ <= 0) return;
  mean_ = mean(value_vec_);
  var_ = var(value_vec_, mean_);
  std::string filepath = base_ + "_statistics.txt";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for SR statistics.");
  ofile << "#Output file for angle calculation with name \'" << name_ << "\'\n";
  ofile << "#Atom group: " << atom_group_->getName() << "\n";
  ofile << "#Average: " << mean_ << " count\n";
  ofile << "#Variance: " << var_ << " count^2\n";
  if(doHistogram){
      std::vector<double> x_vals;
      std::vector<int> y_vals;
      ofile << "#Histogram: nv     count\n";
      makeHistogram(value_vec_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }
  }
  ofile.close();
  if(doTimeseries){
    filepath = base_ + "_timeseries.txt";
    ofile.open(filepath);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for SR timeseries.");
    ofile << "#Output file for nv calculation with name \'" << name_ << "\'\n";   
    ofile << "Timeseries: time (ps)     step     nv\n";
    for(std::size_t i = 0; i < value_vec_.size(); i++){
        ofile << time_vec_[i] << "     " << step_vec_[i] << "     " << value_vec_[i] << "\n"; 
    }
    ofile.close();
  }
  return;
}