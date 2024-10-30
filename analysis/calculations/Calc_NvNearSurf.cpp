#include "Calc_NvNearSurf.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/smearfunctions.hpp"
#include "../../tools/stlmath.hpp"
Calc_NvNearSurf::Calc_NvNearSurf(InputPack& input):Calculation_Histogram{input}
{
  input.params().readString("probe_volume", KeyType::Required, pv_name_);
  auto pv_pointer = input.findProbeVolume(pv_name_);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer;

  std::string agname, surfname, icename;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  
  input.params().readString("surf_group", KeyType::Required, surfname);
  surf_group_ = input.findAtomGroup(surfname);
  FANCY_ASSERT(surf_group_ != 0, "Failed to find specified surf atom group.");

  input.params().readString("ice_group", KeyType::Required, icename);
  ice_group_ = input.findAtomGroup(icename);  
  FANCY_ASSERT(ice_group_ != 0, "Failed to find specified ice atom group.");

  input.params().readNumber("ice_radius", KeyType::Required, ice_radius_);
  input.params().readNumber("ice_threshold_bulk", KeyType::Required, bulk_ice_thresh_);
  input.params().readNumber("ice_threshold_bridging", KeyType::Required, bridging_ice_thresh_);
  input.params().readNumber("ice_smear", KeyType::Required, ice_smear_);
  input.params().readNumber("surface_radius", KeyType::Required, surface_radius_);
  input.params().readNumber("surface_threshold", KeyType::Required, surface_thresh_);   
  input.params().readNumber("surface_smear", KeyType::Required, surface_smear_);

  bOutput = input.params().readString("trajfile", KeyType::Optional, trajfile_);
  if(bOutput){
    trajout_.open(trajfile_);
    FANCY_ASSERT(trajout_.is_open(), "Failed to open output file for trajectory in NvNearSurf");
  }
  return;
}


double Calc_NvNearSurf::neighborEval(const Vec3<double>& pos, CellGrid& cg, double radius, double smear, double cutoff) const{
  double eval=0.0;
  //get neighboring indices
  auto indices = cg.getNearbyIndices(pos);
  //determine how many are within threshold with indus switching function
  for(auto index : indices){
    Vec3<double> pos_temp = box->atoms[index].x;
    //get pbc distance
    double distance = getDistance(pos_temp, pos, box_size_);
    eval += h_r(distance, radius, smear, cutoff);
  }
  return eval;
}

void Calc_NvNearSurf::calculate(){
  if(!doCalculate()) return;
  //put surface atoms and ice-like atoms on neighbor list for computing thresholds
  //c_surf_ and c_ice_

  auto& indices = atom_group_->getIndices();
  auto nidx = indices.size();
  updateNeighborLists();
  ice_neighbors_.clear(); surf_neighbors_.clear(), pv_eval_.clear();
  ice_neighbors_.resize(nidx, 0.0); surf_neighbors_.resize(nidx, 0.0); pv_eval_.resize(nidx, 0.0);
  double sum = 0.0;
  for(int i = 0; i < nidx; i++){
    int idx = atom_group_->getIndices()[i];
    Vec3<double> pos = box->atoms[idx].x; 
    double pv_eval_temp = pv_->compute(pos);
    if(pv_eval_temp  == 0.0) continue; //don't bother if it isn't in the observation volume
    surf_neighbors_[i] = neighborEval(pos, c_surf, surface_radius_, surface_smear_, 2.0*surface_smear_);
    ice_neighbors_[i] = neighborEval(pos, c_ice, ice_radius_, ice_smear_, 2.0*ice_smear_);
    double isf_bridging = 1-h_r(ice_neighbors_[i], bridging_ice_thresh_, 0.01, 0.02); 
    double isf_bulk_bar = h_r(ice_neighbors_[i], bulk_ice_thresh_, 0.01, 0.02); 
    double ssf = 1-h_r(surf_neighbors_[i], surface_thresh_, 0.01, 0.02); 
    //1 if atom doesn't meet surface criteria
    double bridging_eval_bar = 1 - isf_bridging*ssf;
    //0 if atom doesn't meet bulk criteria or surface criteria
    double bulk_eval = 1 - bridging_eval_bar * isf_bulk_bar;
    pv_eval_[i] = pv_eval_temp * bulk_eval;
    sum += pv_eval_[i];
  }

  value_ = sum;

  if(doTimeseries || doHistogram){
    count_vec_.push_back(value_);
    step_vec_.push_back(box->frame_counter);
    time_vec_.push_back(box->time);
  }
  return;
}


std::string Calc_NvNearSurf::printConsoleReport(){
  std::stringstream ss;
  ss << "Name: " << name_ << ", Type: " << "Nv" << ", Output: " << value_ << std::endl;
  return ss.str(); 
}

void Calc_NvNearSurf::finalOutput(){
  if(bOutput) trajout_.close();
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

void Calc_NvNearSurf::printXYZ(){
  auto indices = atom_group_->getIndices();
  auto nidx = indices.size();
  trajout_ << nidx << "\n";
  trajout_ << box->time << "\n";
  for(int i = 0; i < nidx; i++){
    Vec3<double> pos = 10.0*box->atoms[indices[i]].x;
    std::string name = box->atoms[indices[i]].name;
    trajout_ << name << "   " << pos[0] << "   " << pos[1] << "   " << pos[2] << "   " << pv_eval_[i] << "   " << ice_neighbors_[i] << "   " << surf_neighbors_[i] << "\n";
  }
}