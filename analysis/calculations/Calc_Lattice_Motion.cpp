#include "Calc_Lattice_Motion.hpp"
Calc_Lattice_Motion::Calc_Lattice_Motion(InputPack& input):Calculation_Histogram{input}
{
  std::string atom_group_name;
  input.params().readString("atom_group", KeyType::Required, atom_group_name);
  atom_group_ = input.findAtomGroup(atom_group_name);
  FANCY_ASSERT(atom_group_, "Failed to find atom group.");
  //index is a global index of an atom that may or may not be part of the atom group
  input.params().readNumber("index", KeyType::Required, ref_atom_index_);
  ref_atom_index_ -= 1;
  num_nn_ = 4;
  input.params().readNumber("nearest_neighbors", KeyType::Optional, num_nn_); 
  nn_indices_.resize(num_nn_, 0);
  nn_positions_.resize(num_nn_, {0.0});
  nn_distances_.resize(num_nn_, 0.0);
  isInitialized = 0;

}
//overwrites the previous position, returns nearest distance
double Calc_Lattice_Motion::getNearestPeriodicImage(Vec3<double>& pos, const Vec3<double>& ref_pos){
  Vec3<double> dx;
  for(int i = 0; i < 3; i++){
    dx[i] = pos[i] - ref_pos[i];
    if(dx[i] > box_size_[i]/2.0){
      pos[i] -= box_size_[i];
      dx[i] -= box_size_[i];
    } 
    else if(dx[i] < -box_size_[i]/2.0){
      pos[i] += box_size_[i];
      dx[i] += box_size_[i];
    } 
  }
  return std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
}
void Calc_Lattice_Motion::update(){
  Calculation_Histogram::update();
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  if(!doCalculate()) return;
  ref_atom_position_ = box->getAtomPosition(ref_atom_index_);
  if(!isInitialized)
  {
    int nidx = atom_group_->getIndexCount();
    std::vector<double> distances(nidx, std::numeric_limits<double>::max());
    std::vector<Vec3<double> > adjusted_positions(nidx, {0.0});
    for(int i = 0; i < nidx; i++){
      int box_index = atom_group_->getIndex(i);
      if(box_index == ref_atom_index_){
        continue; //don't compare reference with itself
      }
      adjusted_positions[i] = box->getAtomPosition(box_index);
      distances[i] = getNearestPeriodicImage(adjusted_positions[i], ref_atom_position_);
    }

    Vec<double> distance_rankings(num_nn_, std::numeric_limits<double>::max());
    Vec<int> idx_rankings(num_nn_, -1);
    for(int i = 0; i < nidx; i++){
      int box_index = atom_group_->getIndex(i);
      if(box_index == ref_atom_index_){
        continue; //don't compare reference with itself
      }
      for(int j = 0; j < num_nn_; j++){
        if(distances[i] < distance_rankings[j]){
          for(int k = num_nn_-1; k > j; k--){
            distance_rankings[k] = distance_rankings[k-1];
            idx_rankings[k] = idx_rankings[k-1];
          }
          distance_rankings[j] = distances[i];
          idx_rankings[j] = box_index;
          break;
        }
      }
    }
    nn_indices_ = idx_rankings;
    nn_distances_ = distance_rankings;
    for(int i = 0; i < num_nn_; i++){
      nn_positions_[i] = box->getAtomPosition(nn_indices_[i]);
      nn_distances_[i] = getNearestPeriodicImage(nn_positions_[i], ref_atom_position_);
    }
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < num_nn_; j++){
        nn_original_com_[i] += nn_positions_[j][i];
      }
      nn_original_com_[i] *= 1.0/num_nn_;
    }
    isInitialized = 1;
  }
  else{
    for(int i = 0; i < num_nn_; i++){
      nn_positions_[i] = box->getAtomPosition(nn_indices_[i]);
      nn_distances_[i] = getNearestPeriodicImage(nn_positions_[i], ref_atom_position_);
    } 
  }

  return;
}

void Calc_Lattice_Motion::calculate(){
  if(!doCalculate()) return;
  
  Vec3<double> nn_com;
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < num_nn_; j++){
      nn_com[i] += nn_positions_[j][i];
    }
    nn_com[i] *= 1.0/num_nn_;
  }

  Vec3<double> shifted_reference_position = ref_atom_position_;
  for(int i = 0; i < 3; i++){
    shifted_reference_position[i] += nn_original_com_[i] - nn_com[i];
  }

  offsets_.push_back(shifted_reference_position);
  time_vec_.push_back(box->time);
  step_vec_.push_back(box->frame_counter);

  return;
}
std::string Calc_Lattice_Motion::printConsoleReport(){
  return "";
}

void Calc_Lattice_Motion::finalOutput(){
  if(output_freq_ <= 0) return;
  Vec3<double> finalCOM = {0.0};
  //get the final com from time series
  for(int j = 0; j < 3; j++){
    for(int i = 0; i < offsets_.size(); i++){
      finalCOM[j] += offsets_[i][j];
    }
    finalCOM[j] *= 1.0/(double)offsets_.size();
  }
  //shift everything to 0
  for(int i = 0; i < offsets_.size(); i++){
    for(int j = 0; j < 3; j++){
      offsets_[i][j] -= finalCOM[j];
    }
  }
  //calculated offsets from 0
  Vec<double> finalDists(offsets_.size(), 0.0), finalx(offsets_.size(), 0.0), finaly(offsets_.size(), 0.0), finalz(offsets_.size(), 0.0);
  for(int i = 0; i < offsets_.size(); i++){
    for(int j = 0; j < 3; j++){
      finalDists[i] += offsets_[i][j]*offsets_[i][j];
    }
    finalDists[i] = std::sqrt(finalDists[i]);
    finalx[i] += offsets_[i][0];
    finaly[i] += offsets_[i][1];
    finalz[i] += offsets_[i][2];
  }

  double mean = 0.0;
  for(int i = 0; i < finalDists.size(); i++){
    mean += finalDists[i];
  }
  mean *= 1.0/(double)finalDists.size();

  double var = 0.0;
  for(int i = 0; i < finalDists.size(); i++){
    var += (finalDists[i]-mean)*(finalDists[i]-mean);
  }
  var *= 1.0/(double)finalDists.size();

  std::string filepath = base_ + "_statistics.txt";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for Lattice Motion statistics.");
  ofile << "#Output file for Lattice Motion calculation with name \'" << name_ << "\'\n";
  ofile << "#Atom group: " << atom_group_->getName() << "\n";
  ofile << "#Average: " << mean << " count\n";
  ofile << "#Variance: " << var << " count^2\n";
  
  if(doHistogram){
      std::vector<double> x_vals;
      std::vector<int> y_vals;
      ofile << "#Histogram distances: bin     count\n";
      makeHistogram(finalDists, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }

      x_vals.clear();
      y_vals.clear();

      ofile << "#Histogram x: bin     count\n";
      makeHistogram(finalx, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }

      x_vals.clear();
      y_vals.clear();

      ofile << "#Histogram y: bin     count\n";
      makeHistogram(finaly, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }    
  
      x_vals.clear();
      y_vals.clear();

      ofile << "#Histogram z: bin     count\n";
      makeHistogram(finalz, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
      for(std::size_t i = 0; i < x_vals.size(); i++){
        ofile << x_vals[i] << "   " << y_vals[i] << "\n";
      }      
  
  }
  ofile.close();
  if(doTimeseries){
    filepath = base_ + "_timeseries.txt";
    std::ofstream ofile(filepath);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for Lattice Motion timeseries.");
    ofile << "#Output file for Lattice Motion calculation with name \'" << name_ << "\'\n";   
    ofile << "Timeseries: time (ps)     step     nv\n";
    for(std::size_t i = 0; i < finalDists.size(); i++){
        ofile << time_vec_[i] << "     " << step_vec_[i] << "     " << finalDists[i] << "     " << offsets_[i][0] <<  "     " << offsets_[i][1] << "     " << offsets_[i][2] << "\n"; 
    }
    ofile.close();
  }

  return;
}