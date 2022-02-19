#include "Calc_IceID.hpp"
#include "../../tools/cellgrid.hpp"
#include "../../tools/stlmath.hpp"
#include <omp.h>
#include <unordered_set>
Calc_IceID::Calc_IceID(InputPack& input) : Calculation_Histogram{input} {
  filename_ = name_ + ".xtc";
  input.params().readString("filename", KeyType::Optional, filename_);

  input.params().readString("ice_group", KeyType::Required, icename_);
  ice_group_ = input.findAtomGroup(icename_);
  FANCY_ASSERT(ice_group_ != 0, "Failed to find ice atom group.");

  input.params().readString("surface_group", KeyType::Required, surfname_);
  surface_group_ = input.findAtomGroup(surfname_);
  FANCY_ASSERT(surface_group_ != 0, "Failed to find surface atom group.");

  input.params().readString("water_group", KeyType::Required, watername_);
  water_group_ = input.findAtomGroup(watername_);
  FANCY_ASSERT(water_group_ != 0, "Failed to find water atom group.");

  input.params().readString("probe_volume", KeyType::Required, pv_name_);
  auto pv_pointer = input.findProbeVolume(pv_name_);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer;

  input.params().readNumber("ice_radius", KeyType::Required, ice_radius_);
  input.params().readNumber("ice_threshold", KeyType::Required, ice_thresh_);
  input.params().readNumber("surface_radius", KeyType::Required, surface_radius_);
  input.params().readNumber("surface_threshold", KeyType::Required, surf_thresh_);  
  nrep_  = 1;
  input.params().readNumber("iterations", KeyType::Optional, nrep_);  
  initialized_ = 0;
  return;
}

void Calc_IceID::update(){
  if(hasUpdated()) return;
  Calculation::update();
  xdr_time_ = current_time_;
  xdr_step_ = current_frame_;
  xdr_natoms_ = box->atoms.size();
  if(!initialized_){
    xdr_prec_ = 1000;
    xdr_x_ = new xdr::rvec[xdr_natoms_];
    output_handle_ = xdr::xdrfile_open(filename_.c_str(), "w");
    FANCY_ASSERT(output_handle_!=0, "xtc output file is null");

    index_out_.open(name_ + ".index");
    initialized_ = 1;
  }
  return;
}

void Calc_IceID::performIteration(std::vector<int>& ice_indices){
	std::unordered_set<int> final_ice_set;
  for(auto index : ice_indices){
    final_ice_set.insert(index);
  }
 //create cellgrid
  Vec3<double> box_vec = {box->boxvec[0][0], box->boxvec[1][1], box->boxvec[2][2]};
	CellGrid c_ice(ice_radius_, box_vec);
  CellGrid c_surf(surface_radius_, box_vec);
  auto& surf_indices = surface_group_->getIndices();
  auto& water_indices = water_group_->getIndices();
  auto& atoms = box->atoms;

  for(int i = 0; i < ice_indices.size(); i++){
    c_ice.addIndexToGrid(ice_indices[i], atoms[ice_indices[i]].x);
  }
  for(int i = 0; i < surf_indices.size(); i++){
    c_surf.addIndexToGrid(surf_indices[i], atoms[surf_indices[i]].x);
  }
  std::vector<int> ice_neighbor_counts(water_indices.size(), 0);
  std::vector<int> surface_neighbor_counts(water_indices.size(), 0);
  //count the number of neighbors of each atom type
  #pragma omp parallel for
  for(int i = 0; i < water_indices.size(); i++){
    int w_idx = water_indices[i];
    if(final_ice_set.find(w_idx) != final_ice_set.end()) continue; //if it's already icelike, don't bother
    auto ice_neighbors = c_ice.getNearbyIndices(water_indices[i], atoms[water_indices[i]].x);
    auto surf_neighbors = c_surf.getNearbyIndices(water_indices[i], atoms[water_indices[i]].x);
    if(ice_neighbors.size() > 0){
      for(auto it = ice_neighbors.begin(); it != ice_neighbors.end();){
        int i_idx = *it;
        double dist = getDistance(atoms[i_idx].x, atoms[w_idx].x, box_vec);
        if(dist > ice_radius_){
          it = ice_neighbors.erase(it);
        }
        else{
          ++it;
        }
      }
    }
    if(surf_neighbors.size() > 0){
      for(auto it = surf_neighbors.begin(); it != surf_neighbors.end();){
        int s_idx = *it;
        double dist = getDistance(atoms[s_idx].x, atoms[w_idx].x, box_vec);
        if(dist > surface_radius_){
          it = surf_neighbors.erase(it);
        }
        else{
          ++it;
        }
      }
    }
    int numIceNeighbors = ice_neighbors.size();
    int numSurfNeighbors = surf_neighbors.size();
    ice_neighbor_counts[i] = numIceNeighbors;
    surface_neighbor_counts[i] = numSurfNeighbors;
  }

  //now for each water index, we have then number of surface and ice neighbors
  for(int i = 0; i < water_indices.size(); i++){
    if(final_ice_set.find(water_indices[i]) != final_ice_set.end()) continue;
    if(ice_neighbor_counts[i] >= ice_thresh_ && surface_neighbor_counts[i] >= surf_thresh_){
      ice_indices.push_back(water_indices[i]);
    }
  }
  return;
}

void Calc_IceID::calculate(){
  if(!doCalculate()) return;
  auto& atoms = box->atoms;
  auto& water_indices = water_group_->getIndices();
  std::unordered_set<int> wi(water_indices.begin(), water_indices.end());
  final_ice_indices_ = ice_group_->getIndices();
  n_original_ = 0.0;
  for(auto index : final_ice_indices_){
    if(wi.find(index) != wi.end()) n_original_ += pv_->compute(atoms[index].x); //atom must be a water to count towards the final tally
  }

  if(nrep_ > 0){
    for(int i = 0; i < nrep_; i++){
      performIteration(final_ice_indices_);
    }
  }

  else{
    int counter = 0;
    int lastcount = final_ice_indices_.size();
    while(counter < 10){
      performIteration(final_ice_indices_);
      if(lastcount == final_ice_indices_.size()) break;
      counter++;
    }
  }

  int n_ = 0;
  for(auto index : final_ice_indices_){
    if(wi.find(index) != wi.end()) n_ += pv_->compute(atoms[index].x); //atom must be a water to count towards the final tally
  }
  t_vec_.push_back(current_time_);
  step_vec_.push_back(current_frame_);
  n_vec_.push_back(n_);
  no_vec_.push_back(n_original_);
  return;
}

void Calc_IceID::output(){
  if(!doOutput()) return;
  for(int i = 0; i < box->atoms.size(); i++){
    for(int j = 0; j < 3; j++){
      xdr_x_[i][j] = -1.0f;
    }
  }

  auto indices = final_ice_indices_;
  for(int i = 0; i < indices.size(); i++){
    int index = indices[i];
    for(int j = 0; j < 3; j++){
      xdr_x_[index][j] = box->atoms[index].x[j];
    }
  }

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
        xdr_box_[i][j] = box->boxvec[i][j];
    }
  }

  xdr::write_xtc(output_handle_, xdr_natoms_, xdr_step_, xdr_time_, xdr_box_, xdr_x_, xdr_prec_);

  index_out_ << current_time_ << "   ";
  for(auto index : final_ice_indices_){
    index_out_ << index << "   ";
  }
  index_out_ << "\n";

  return;
}

void Calc_IceID::finalOutput(){
  
  double n_avg = mean(n_vec_);
  double n_var = var(n_vec_, n_avg);
  double no_avg = mean(no_vec_);
  double no_var = var(no_vec_, no_avg);
  std::string filepath = base_ + "_statistics.txt";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for IceID");
  ofile << "#Output file for IceID calculation with name \'" << name_ << "\'\n";
  ofile << "#Atom group: " << water_group_->getName() << "\n";
  ofile << "#Average: " << n_avg << " count, Original Average: " << no_avg << "\n";
  ofile << "#Variance: " << n_var << " count^2, Original Variance: " << no_var << "\n";
  if(doHistogram){
      std::vector<double> x_vals;
      std::vector<int> y_vals;
      ofile << "#Histogram: nv     count\n";
      makeHistogram(n_vec_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }
  }
  ofile.close();
  if(doTimeseries){
    filepath = base_ + "_timeseries.txt";
    std::ofstream ofile(filepath);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for nv timeseries.");
    ofile << "#Output file for IceID calculation with name \'" << name_ << "\'\n";   
    ofile << "Timeseries: time (ps)     step     nv     original count\n";
    for(std::size_t i = 0; i < n_vec_.size(); i++){
        ofile << t_vec_[i] << "     " << step_vec_[i] << "     " << n_vec_[i] << "\n"; 
    }
    ofile.close();
  }

  xdr::xdrfile_close(output_handle_);
  delete xdr_x_;
  index_out_.close();
  return;
}