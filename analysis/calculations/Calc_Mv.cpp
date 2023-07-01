#include "Calc_Mv.hpp"
Calc_Mv::Calc_Mv(InputPack& input):Calculation_Histogram{input}
{
  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");

  double qbar_max, qbar_cut, qbar_sigma;
  int qlmi_l;
  double qlmi_rmax, qlmi_cut, qlmi_sigma;
  double nn_min, nn_max, nn_cut, nn_sigma;
  Vec3<double> pv_min, pv_max;
  double pv_sigma, pv_cut;
  Vec<double> xr, yr, zr;
  input.params().readVector("x_range", KeyType::Required, xr);
  input.params().readVector("y_range", KeyType::Required, yr);
  input.params().readVector("z_range", KeyType::Required, zr);
  input.params().readNumber("pv_sigma", KeyType::Required, pv_sigma);
  input.params().readNumber("pv_cut", KeyType::Required, pv_cut);  
  pv_min[0] = xr[0]; pv_max[0] = xr[1];
  pv_min[1] = yr[0]; pv_max[1] = yr[1];
  pv_min[2] = zr[0]; pv_max[2] = zr[1];

  input.params().readNumber("qbar_cutoff", KeyType::Required, qbar_cut);
  input.params().readNumber("qbar_sigma", KeyType::Required, qbar_sigma);
  input.params().readNumber("qbar_max", KeyType::Required, qbar_max);  

  input.params().readNumber("qlmi_order", KeyType::Required, qlmi_l);
  input.params().readNumber("qlmi_cutoff", KeyType::Required, qlmi_cut);
  input.params().readNumber("qlmi_sigma", KeyType::Required, qlmi_sigma);
  input.params().readNumber("qlmi_rmax", KeyType::Required, qlmi_rmax);   
  qlmi_rmax_ = qlmi_rmax;
  qlmi_rcut_ = qlmi_cut;
  qlmi_order_ = qlmi_l;
  input.params().readNumber("nn_min", KeyType::Required, nn_min);
  input.params().readNumber("nn_max", KeyType::Required, nn_max);
  input.params().readNumber("nn_cut", KeyType::Required, nn_cut);
  input.params().readNumber("nn_sigma", KeyType::Required, nn_sigma);
  //
  input.params().readNumber("dump", KeyType::Required, dump_frame_);
  //simpleTwoSidedSwitch* h4nn_switch_;
  //simpleSwitch* qbar_cutoff_;
  //cuboidalINDUSVolume* probe_volume_;
  qlmi_ = new qlmi(qlmi_l, qlmi_rmax, qlmi_sigma, qlmi_cut);
  h4nn_switch_ = new simpleTwoSidedSwitch(nn_min, nn_max, nn_sigma, nn_cut);
  qbar_cutoff_ = new simpleSwitch(qbar_max, qbar_sigma, qbar_cut);
  probe_volume_ = new cuboidalINDUSVolume(pv_min, pv_max, pv_sigma, pv_cut);
  Vec3<double> shell_min = pv_min;
  Vec3<double> shell_max = pv_max;  
  for(int i = 0; i < 3; i++){
    shell_min[i] -= 2.0*(qlmi_rmax + qlmi_cut);
    shell_max[i] += 2.0*(qlmi_rmax + qlmi_cut);
  }
  shell_ = new cuboidalINDUSVolume(shell_min, shell_max, pv_sigma, pv_cut);
  input.params().readNumber("dump", KeyType::Optional, dump_frame_);

  //for(double i = 0; i < 0.5; i+=0.01){
  //  Vec3<double> pos = {i,0,0};
  //  std::cout << i << "  " << qlmi_->qlmi_test_probe(pos, {0,0,0}) << "\n";
  //}
  return;
}
void Calc_Mv::calculate(){
  if(!doCalculate()) return;
  //update all box_volumes
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  //std::cout << "legendre 1,2 " << p_legendre(2,1,0.3) << std::endl;
  probe_volume_->update_box_size(box_size_);
  shell_->update_box_size(box_size_);
  qlmi_->updateBoxSize(box_size_);
  double sum = 0.0;
  int natoms = atom_group_->getIndices().size();
  Vec<Vec3<double> > positions(natoms);
  for(int i = 0; i < natoms; i++){
    int idx = atom_group_->getIndices()[i];
    positions[i] = box->atoms[idx].x;
  }
  //for each atom get htilde_shell, if htilde_shell > 0, include positions in final positions vector
  Vec<Vec3<double>> final_positions;
  final_positions.reserve(natoms);
  for(const auto& position : positions){
    double h = shell_->calculate(position);
    if(h > 0) final_positions.push_back(position);
  }
  int nfinal = final_positions.size();
  //build neighbor list with qlmi
  buildCellGrid(final_positions);
  //loop through all particles, get their neighbors, and compute qlmi
  Vec<Vec<std::complex<double> > > qlm_mat(final_positions.size());
  Vec<double> nn_vec(final_positions.size());
  Vec<Vec<int> > ni_vec(final_positions.size(), Vec<int>());
  Vec<Vec<double> > htilde_neighbors_vec(final_positions.size(),  Vec<double>());
  for(int i = 0; i < final_positions.size(); i++){
    Vec<Vec3<double> > neighbor_positions;
    Vec<int> neighbor_indices, neighbor_indices2;
    getNeighbors(i, final_positions[i], neighbor_indices);
    for(auto index : neighbor_indices){
      double distance = getDistance(final_positions[index], final_positions[i], box_size_);
      if(distance < qlmi_rmax_ + qlmi_rcut_){
        neighbor_positions.push_back(final_positions[index]);
        neighbor_indices2.push_back(index);
      }
    }
    ni_vec[i] = neighbor_indices2;
    std::vector<std::complex<double> > ylm_sum_vec;
    Vec<double> htilde_neighbors;
    double nn_i = qlmi_->compute(final_positions[i], neighbor_positions, ylm_sum_vec, htilde_neighbors);
    htilde_neighbors_vec[i] = htilde_neighbors;
    qlm_mat[i] = ylm_sum_vec;
    nn_vec[i] = nn_i;
  }
  //now have qlm(i) values for all indices that could potentially be relevant
  //just need to average over neighbors and use switching functions to get qbar values
  std::vector<double> hi_values(final_positions.size());
  std::vector<double> hnn_values(final_positions.size());
  std::vector<double> qbar_values(final_positions.size());
  std::vector<double> mtilde_values(final_positions.size());
  for(int i = 0; i < final_positions.size(); i++){
    hi_values[i] = probe_volume_->calculate(final_positions[i]);
    hnn_values[i] = h4nn_switch_->calculate(nn_vec[i]);
    if(hi_values[i] > 0.0){
      double m_sum = 0.0;
      for(int m = 0; m < 2*qlmi_order_ + 1; m++){
        std::complex<double> hviqlm = 0.0;
        for(int k = 0; k < htilde_neighbors_vec[i].size(); k++){
          hviqlm += qlm_mat[ni_vec[i][k]][m]*htilde_neighbors_vec[i][k];
        }
        std::complex<double> qbarmi = (1.0/(1.0+nn_vec[i])) * (qlm_mat[i][m] + hviqlm);
        m_sum += std::pow(std::abs(qbarmi),2);
      }
      qbar_values[i] = sqrt(m_sum * (4.0*M_PI/13.0));
      mtilde_values[i] = qbar_cutoff_->calculate(qbar_values[i]) * hnn_values[i] * hi_values[i];
    }
    else{
      qbar_values[i] = 0.0;
      mtilde_values[i] = 0.0;
    }
  }
  sum = 0.0;
  for(int i = 0; i < mtilde_values.size(); i++){
    sum += mtilde_values[i];
  }
  value_ = sum;

  if(atomwiseNeeded){
    positions_ = final_positions;
    mvals_ = mtilde_values;
  }

  if(doTimeseries || doHistogram){
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  if(current_frame_ == dump_frame_){
    std::ofstream ofile("mv_dump.txt");
    ofile << "#DUMP FILE FOR CALC_MV" << std::endl;
    for(int i = 0; i < final_positions.size(); i++){
      ofile << i << "  " << hi_values[i] <<  "  " << hnn_values[i] <<  "  " << qbar_values[i] << "  " << mtilde_values[i] << "\n";
    }
    ofile.close();
  }
  return;
}
std::string Calc_Mv::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Nv" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}

void Calc_Mv::finalOutput(){
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
