#include "Calc_Angle.hpp"
Calc_Angle::Calc_Angle(InputPack& input):Calculation{input}
{
  Vec3<std::string> agnames;
  input.params().readString("group1", KeyType::Required, agnames[0]);
  input.params().readString("group2", KeyType::Required, agnames[1]);
  input.params().readString("group3", KeyType::Required, agnames[2]);
  FANCY_ASSERT(agnames[0] != agnames[1] && agnames[1] != agnames[2] && agnames[0] != agnames[3],
      "Identical atomgroups given for two parts of angle calculation.");
  for(int i = 0; i < 3; i++){
      atom_groups_[i] = input.findAtomGroup(agnames[i]);
      FANCY_ASSERT(atom_groups_[i] != 0, "Failed to find atom group.");
  }
  FANCY_ASSERT( (atom_groups_[0]->getIndices().size() == atom_groups_[1]->getIndices().size()) && 
                (atom_groups_[0]->getIndices().size() == atom_groups_[2]->getIndices().size()), 
                "Mismatched atom group sizes.");
  return;
}
inline double f_angle(Atom a1, Atom a2, Atom a3){
    Vec3<double> v21,v23;
    double num = 0, norm1 = 0, norm2 = 0;
    for(int i = 0; i < 3; i++){
        v21[i] = a1.x[i] - a2.x[i];
        v23[i] = a3.x[i] - a2.x[i];
        num += v21[i]*v23[i];
        norm1 += v21[i]*v21[i];
        norm2 += v23[i]*v23[i];
    }
    double denom = sqrt(norm1)*sqrt(norm2);
    return acos(num/denom);
}
void Calc_Angle::calculate(){
  if(!doCalculate()) return;
  int natoms = atom_groups_[0]->getIndices().size();
  Vec<double> angle_distribution_step(natoms);
  value_ = 0.0;
  for(int i = 0; i < natoms; i++){
    Vec3<int> idx;
    for(int j = 0; j<3; j++){
        idx[j] = atom_groups_[j]->getIndices()[i];
    }
    Atom a1 = box->atoms[idx[0]];
    Atom a2 = box->atoms[idx[0]];
    Atom a3 = box->atoms[idx[0]];
    double angle = f_angle(a1, a2, a3);
    value_ += angle;
    angle_distribution_step[i] = angle;
  }
  value_ *= 1.0/(double)natoms;
  angle_distribution_.push_back(angle_distribution_step);
  if(doTimeseries){
    angle_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  //this calculation has no per-timestep output, so print output function is blank (defined in Calculations.hpp)
  if(doOutput()) printOutput();
  return;
}
std::string Calc_Angle::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "angle" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}

void Calc_Angle::finalOutput(){
  if(output_freq_ <= 0) return;
  //get average angle by looping through distribution
  //should be conducive to an omp reduction, so not optimizing up front
  double sum = 0.0;
  int counts = 0;
  for(std::size_t i = 0; i < angle_distribution_.size(); i++){
      for(std::size_t j = 0; j < angle_distribution_[i].size(); j++){
        sum += angle_distribution_[i][j];
      }
      counts += angle_distribution_[i].size();
  }
  mean_ = sum/(double)counts;

  double sum2 = 0;
  for(std::size_t i = 0; i < angle_distribution_.size(); i++){
      for(std::size_t j = 0; j < angle_distribution_[i].size(); j++){
        sum2 += pow(angle_distribution_[i][j]-mean_, 2);
      }
  }
  var_ = sum2/(double)counts;
  std::string filepath = base_ + "_statistics.txt";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for angle statistics.");
  ofile << "#Output file for angle calculation with name \'" << name_ << "\'\n";
  ofile << "#Atom groups: " << atom_groups_[0]->getName() << " " << atom_groups_[1]->getName() << " " << atom_groups_[2]->getName() << "\n";
  ofile << "#Average: " << mean_ << " radians\n";
  ofile << "#Variance: " << var_ << " radians^2\n";
  if(doHistogram){
      std::vector<double> x_vals;
      std::vector<int> y_vals;
      std::cout << "#Histogram: angle (rads)     count\n";
      makeHistogram2d(angle_distribution_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }
  }
  ofile.close();
  if(doTimeseries){
    filepath = base_ + "_timeseries.txt";
    std::ofstream ofile(filepath);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for angles timeseries.");
    std::cout << "#Output file for angle calculation wityh name \'" << name_ << "\'\n";   
    std::cout << "Timeseries: time (ps)     step     angle (rads)\n";
    for(std::size_t i = 0; i < angle_vec_.size(); i++){
        std::cout << time_vec_[i] << "     " << step_vec_[i] << "     " << angle_vec_[i] << "\n"; 
    }
    ofile.close();
  }
  return;
}