#include "Calc_HPolarity.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/pbcfunctions.hpp"
Calc_HPolarity::Calc_HPolarity(InputPack& input):Calculation_Histogram{input}
{
  input.params().readString("probe_volume", KeyType::Required, pv_name_);
  auto pv_pointer = input.findProbeVolume(pv_name_);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer;
  //right now this is only set up to work for systems where two atoms exist along
  //the dipole, for tip4p this would be the OW and MW atoms
  std::string h_agname, t_agname;
  input.params().readString("head_atom_group", KeyType::Required, h_agname);
  input.params().readString("tail_atom_group", KeyType::Required, t_agname);
  t_group_ = input.findAtomGroup(t_agname);
  FANCY_ASSERT(t_group_ != 0, "Failed to find tail atom group.");
  h_group_ = input.findAtomGroup(h_agname);
  FANCY_ASSERT(h_group_ != 0, "Failed to find head atom group.");

  std::vector<double> surface_normal;
  input.params().readVector("surface_normal", KeyType::Required, surface_normal);
  FANCY_ASSERT(surface_normal.size() == 3, "surface normal must be a 3d vector");
  //normalize to 1
  surface_normal = surface_normal * (1.0/norm2(surface_normal));
  for(int i = 0; i < 3; i++){
    normal_[i] = surface_normal[i];
  }
  return;
}
void Calc_HPolarity::calculate(){
  if(!doCalculate()) return;
  float sum = 0.0;
  double pdotn = 0.0, abspdotn=0.0;
  double nv=0.0;
  #pragma omp parallel for reduction(+:pdotn) reduction(+:abspdotn) reduction(+:nv)
  for(int i = 0; i < h_group_->getIndices().size(); i++){
    int idx = h_group_->getIndices()[i];
    int idx2 = t_group_->getIndices()[i];
    Vec3<double> x = box->atoms[idx].x;
    Vec3<double> xt = box->atoms[idx2].x;
    double htemp = pv_->compute(box->atoms[idx].x);
    if(htemp > 0.5){
      getNearestImage3D(xt, x, box_size_);
      double pdotn_temp = dot((x - xt), normal_);
      pdotn += pdotn_temp;
      abspdotn += fabs(pdotn_temp);
      nv += htemp;
    }
  }
  pdotn /= abspdotn;
  value_ = pdotn;
  if(doTimeseries || doHistogram){
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  return;
}
std::string Calc_HPolarity::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Xi" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}

void Calc_HPolarity::finalOutput(){
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
  ofile << "#Head atom group: " << h_group_->getName() << "\n";
  ofile << "#Tail atom group: " << t_group_->getName() << "\n";
  ofile << "#Surface normal: " << "[ " << normal_[0] << "  " << normal_[1] << "  " << normal_[2] << " ]" << "\n";
  ofile << "#Average: " << mean_ << " count\n";
  ofile << "#Variance: " << var_ << " count^2\n";
  if(doHistogram){
      std::vector<double> x_vals;
      std::vector<int> y_vals;
      ofile << "#Histogram: xi     count\n";
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
    ofile << "#Output file for xi calculation with name \'" << name_ << "\'\n";   
    ofile << "Timeseries: time (ps)     step     nv\n";
    for(std::size_t i = 0; i < count_vec_.size(); i++){
        ofile << time_vec_[i] << "     " << step_vec_[i] << "     " << count_vec_[i] << "\n"; 
    }
    ofile.close();
  }
  return;
}
