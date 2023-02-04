#include "Calc_Nv_SWIPES.hpp"
Calc_Nv_SWIPES::Calc_Nv_SWIPES(InputPack& input):Calculation_Histogram{input}
{
  auto& params = input.params();
  std::string agname;
  params.readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group."); 
  //in order, [ xmin dx xmax ]
  int dim = 0;
  Vec<double> template_min, template_max;
  params.readNumber("dim", ParameterPack::KeyType::Required, dim);
  params.readVector("template_min", ParameterPack::KeyType::Required, template_min);
  params.readVector("template_max", ParameterPack::KeyType::Required, template_max);
  FANCY_ASSERT(dim >= 0 && dim < 3, "dim should correspond to the axis that you're incrementing over and should be between 0 and 2");
  Vec<double> xmin_range;
  Vec<double> xmax_range;
  params.readVector("xmin", ParameterPack::KeyType::Required, xmin_range);
  params.readVector("xmax", ParameterPack::KeyType::Required, xmax_range);
  FANCY_ASSERT(xmin_range.size() == 3, "Invalid size given for xmin_range, you provided " 
              + std::to_string(xmin_range.size()) + " expected 3, [ min increment max ]");
  FANCY_ASSERT(xmax_range.size() == 3, "Invalid size given for xmax_range, you provided " 
              + std::to_string(xmax_range.size()) + " expected 3, [ min increment max ]");
  //should form an MINxMAX matrix of entries with every combination of origin and spacing
  FANCY_ASSERT(xmin_range[0] < xmin_range[2], "Invalid range provided, for arguments [ min increment max ] \\
  xmin should be less than xmax");
  FANCY_ASSERT(xmax_range[0] < xmax_range[2], "Invalid range provided, for arguments [ min increment max ] \\
  xmin should be less than xmax");  
  FANCY_ASSERT(xmin_range[1] > 0, "increment should be > 0 for xmin");
  FANCY_ASSERT(xmax_range[1] > 0, "increment should be > 0 for xmax");

  for(double x1 = xmin_range[0]; x1 <= xmin_range[2]; x1+=xmin_range[1]){
    for(double x2 = xmax_range[0]; x2 <= xmax_range[2]; x2+=xmax_range[1]){
      std::array<double,2> xrange = {template_min[0], template_max[0]}, 
                           yrange = {template_min[1], template_max[1]},
                           zrange = {template_min[2], template_max[2]};
      if(dim == 0) xrange = {x1, x2}; if(dim == 1) yrange = {x1, x2}; if(dim == 2) zrange = {x1, x2};
      PV_DiscreteRect pv_temp = PV_DiscreteRect(xrange, yrange, zrange);
      pv_set_.push_back(pv_temp);
      ranges_.push_back({x1,x2});
    }
  }
  //should have a range of pv sizes now, just have to store the 
  return;
}
void Calc_Nv_SWIPES::calculate(){
  if(!doCalculate()) return;
  double sum = 0.0;
  std::vector<double> data_step(pv_set_.size());
  for(int i = 0; i < atom_group_->getIndices().size(); i++){
    int idx = atom_group_->getIndices()[i];
    auto& x = box->atoms[idx].x;
    for(int j = 0; j < pv_set_.size(); j++){
      data_step[j] = pv_set_[j].compute(x);
      sum += data_step[j] / (ranges_[j][1]-ranges_[j][0]);
    }
  }
  value_ = sum / (double)data_step.size();
  if(doTimeseries){
    ts_values_.push_back(data_step);
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  return;
}
std::string Calc_Nv_SWIPES::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Nv" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}

void Calc_Nv_SWIPES::finalOutput(){
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
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for nv statistics.");
  ofile << "#Output file for nv calculation with name \'" << name_ << "\'\n";
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
    filepath = base_ + "_nvSWIPES.txt";
    ofile.open(filepath);
    std::vector<double> averages(pv_set_.size(),0.0);
    for(std::size_t i = 0; i < ts_values_.size(); i++){
      ofile << time_vec_[i] << "  " << step_vec_[i] << "  ";
      for(int j = 0; j < ts_values_[i].size(); j++){
        ofile << ts_values_[i][j] << "  ";;
        averages[j] += ts_values_[i][j];
      }
      ofile << "\n";
    }
    ofile << "#AVERAGES ";
    for(auto& val : averages){
      val *= 1.0/pv_set_.size();
      ofile << val << "  ";
    }
    ofile.close();
  }
  return;
}
