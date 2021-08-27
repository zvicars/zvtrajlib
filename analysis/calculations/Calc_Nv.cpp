#include "Calc_Nv.hpp"
Calc_Nv::Calc_Nv(InputPack& input):Calculation{input}
{
  input.params().readString("probe_volume", KeyType::Required, pv_name_);
  auto pv_pointer = input.findProbeVolume(pv_name_);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer;

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");

  input.params().readNumber("dump", KeyType::Optional, dump_frame_);

  return;
}
void Calc_Nv::calculate(){
  if(!doCalculate()) return;
  float sum = 0.0;
  for(int i = 0; i < atom_group_->getIndices().size(); i++){
    int idx = atom_group_->getIndices()[i];
    sum += pv_->compute(box->atoms[idx].x);
  }
  value_ = sum;
  if(doTimeseries || doHistogram){
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  return;
}
std::string Calc_Nv::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Nv" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}

void Calc_Nv::finalOutput(){
  if(output_freq_ <= 0) return;
  double sum = 0.0;
  int counts = 0;
  for(std::size_t i = 0; i < count_vec_.size(); i++){
    sum += count_vec_[i];
    counts++;
  }
  mean_ = sum/(double)counts;

  double sum2 = 0.0;
  for(std::size_t i = 0; i < count_vec_.size(); i++){
    sum2 += pow(count_vec_[i]-mean_, 2);
  }
  var_ = sum2/(double)counts;
  std::string filepath = base_ + "_statistics.txt";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for angle statistics.");
  ofile << "#Output file for angle calculation with name \'" << name_ << "\'\n";
  ofile << "#Atom group: " << atom_group_->getName() << "\n";
  ofile << "#Average: " << mean_ << " count\n";
  ofile << "#Variance: " << var_ << " count^2\n";
  if(doHistogram){
      std::vector<double> x_vals;
      std::vector<int> y_vals;
      ofile << "#Histogram: nv     count\n";
      makeHistogram(count_vec_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }
  }
  ofile.close();
  if(doTimeseries){
    filepath = base_ + "_timeseries.txt";
    std::ofstream ofile(filepath);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for nv timeseries.");
    ofile << "#Output file for nv calculation with name \'" << name_ << "\'\n";   
    ofile << "Timeseries: time (ps)     step     nv\n";
    for(std::size_t i = 0; i < count_vec_.size(); i++){
        ofile << time_vec_[i] << "     " << step_vec_[i] << "     " << count_vec_[i] << "\n"; 
    }
    ofile.close();
  }
  return;
}
