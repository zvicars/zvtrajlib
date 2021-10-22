#include "Calc_Dipole.hpp"
#include "../helper/make_histogram.hpp"
Calc_Dipole::Calc_Dipole(InputPack& input):Calculation{input}
{
  std::string agname;
  input.params().readString("atomgroup", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  Vec<double> normal;
  input.params().readVector("normal", KeyType::Required, normal); //surface normal for cos theta calc
  FANCY_ASSERT(normal.size() == 3, "Improperly sized normal vector in dipole calc");
  double norm = 0.0;
  for(int i = 0; i < 3; i++){
    norm += normal[i]*normal[i];
  }
  norm = sqrt(norm);
  for(int i = 0; i < 3; i++){
    normal_[i] = normal[i]/norm;
  }
  FANCY_ASSERT(atom_group_ != 0, "Failed to find atom group.");
  FANCY_ASSERT(atom_group_->getType() == "resname", "Can only use resname type at the moment.");
  auto molinfo = input.params().findParameterPacks("molecule_info", KeyType::Required);
  FANCY_ASSERT(molinfo.size() == 1, "Expecting only a single molinfo object in the calculation");
  molinfo[0]->readNumber("natoms", KeyType::Required, at_per_mol_);
  molinfo[0]->readVector("charges", KeyType::Required, charges_);
  molinfo[0]->readVector("masses", KeyType::Required,  masses_);
  molinfo[0]->readVector("include", KeyType::Required, included_atoms_);
  FANCY_ASSERT(charges_.size() == at_per_mol_, "Improperly sized charge vector");
  FANCY_ASSERT(masses_.size() == at_per_mol_, "Improperly sized mass vector");
  FANCY_ASSERT(included_atoms_.size() == at_per_mol_, "Improperly sized included atoms");
  costheta_avg_total_counter_ = 0;
  costheta_avg_total_ = 0.0; 
  return;
}


void Calc_Dipole::calculate(){
  auto indices = atom_group_->getIndices();
  int mol_it = 0, at_it = 0, mol_at_it=0;
  Vec3<double> dipole_vec = {0.0};
  while(mol_it < nmols_)
  {
    auto pos = box->atoms[indices[at_it]].x;
    if(included_atoms_[mol_at_it]){
        for(int i = 0; i < 3; i++){
            dipole_vec[i] += charges_[i]*pos[i];
        }
    }
    if(at_it%at_per_mol_ == at_per_mol_-1){
      dipoles_[mol_it] = dipole_vec;
      mol_it++;
      mol_at_it = 0;
    }
    at_it++;
    mol_at_it++;
  }
  //should have calculated an axis for each molecule, now calculate the cos(theta) between that and the surface normal
  for(int i = 0; i < dipoles_.size(); i++){
    double norm = 0.0;
    for(int j = 0; j < 3; j++){
      norm += dipoles_[i][j]*dipoles_[i][j];
    }
    double invnorm = 1.0/std::sqrt(norm);
    double costheta = invnorm*(dipoles_[i][0]*normal_[0] + dipoles_[i][1]*normal_[1] + dipoles_[i][2]*normal_[2]);
    costhetas_[i] = costheta;
    costheta_avg_frame_ += costheta;
    costheta_avg_total_ += costheta;
    costheta_avg_total_counter_++;
  }
    costheta_avg_frame_ *= 1.0/(double)nmols_;
    costhetas_.push_back(costheta_avg_frame_);
    times_.push_back(box->time);
    steps_.push_back(box->frame);
  return;
}

void Calc_Dipole::update(){
  Calculation::update();
  auto indices = atom_group_->getIndices();
  natoms_ = indices.size();
  nmols_ = natoms_/at_per_mol_;
  dipoles_.resize(nmols_);
  costhetas_.resize(nmols_);
  costheta_avg_frame_ = 0.0;
  return; 
}
std::string Calc_Dipole::printConsoleReport(){
  return "";
}

void Calc_Dipole::finalOutput(){
  if(output_freq_ <= 0) return;
  double sum = 0.0;
  int counts = 0;
  for(std::size_t i = 0; i < costhetas_.size(); i++){
    sum += costhetas_[i];
    counts++;
  }
  costheta_avg_total_ *= 1.0/(double)costheta_avg_total_counter_;

  double sum2 = 0.0;
  for(std::size_t i = 0; i < costhetas_.size(); i++){
    sum2 += pow(costhetas_[i]-costheta_avg_total_, 2);
  }
  var_ = sum2/(double)counts;
  std::string filepath = base_ + "_statistics.txt";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for angle statistics.");
  ofile << "#Output file for angle calculation with name \'" << name_ << "\'\n";
  ofile << "#Atom group: " << atom_group_->getName() << "\n";
  ofile << "#Average: " << costheta_avg_total_ << " count\n";
  ofile << "#Variance: " << var_ << " count^2\n";
  if(doHistogram){
      std::vector<double> x_vals;
      std::vector<int> y_vals;
      ofile << "#Histogram: nv     count\n";
      makeHistogram(costhetas_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
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
    for(std::size_t i = 0; i < costhetas_.size(); i++){
        ofile << times_[i] << "     " << steps_[i] << "     " << costhetas_[i] << "\n"; 
    }
    ofile.close();
  }
  return;
}