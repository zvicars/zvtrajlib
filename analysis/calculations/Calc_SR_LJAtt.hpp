//Zachariah Vicars 7/7/2022
//Calc SR LJAtt computes the attractive part of the lennard jones potential

#pragma once
#include "Calc_SR.hpp"
#include "../../tools/SimpleForcefields.hpp"
#include "../../tools/Matrix.hpp"
#include <set>
#include <limits>
class Calc_SR_LJAtt : public Calc_SR{
public:
  Calc_SR_LJAtt(InputPack& input):Calc_SR{input}{
    input.params().readVector("sigmas", KeyType::Required, sigmas_);
    input.params().readVector("epsilons", KeyType::Required, epsilons_);
    //pairs of atoms specified with corresponding entries in vector
    input.params().readVector("atom_names1", KeyType::Required, atom_names_);
    input.params().readVector("atom_names2", KeyType::Required, atom_names2_);
    FANCY_ASSERT(atom_names_.size() == sigmas_.size(), "Need one sigma per atom name.");
    FANCY_ASSERT(atom_names_.size() == epsilons_.size(), "Need one epsilon per atom name.");
    FANCY_ASSERT(atom_names_.size() == atom_names2_.size(), "Need one epsilon per atom name.");
    std::set<std::string> all_names;
    for(auto name : atom_names_){
      all_names.insert(atom_names_.begin(), atom_names_.end()); 
      all_names.insert(atom_names2_.begin(), atom_names2_.end()); 
    }
    std::size_t total_unique_names=0;
    for(auto entry : all_names){
      index_map_.insert({entry, total_unique_names});
      total_unique_names++;
    }
    sigma16_prefactor_ = 1.122462048;
    epsilon_matrix_.initialize({total_unique_names, total_unique_names});
    sigma_matrix_.initialize({total_unique_names, total_unique_names});
    epsilon_matrix_.fill(0.0);
    sigma_matrix_.fill(0.0);
    
    for(int i = 0; i < atom_names_.size(); i++){
      auto name1 = atom_names_[i];
      auto name2 = atom_names2_[i];
      std::size_t nidx1 = index_map_.find(name1)->second;
      std::size_t nidx2 = index_map_.find(name2)->second;
      epsilon_matrix_.at((std::array<std::size_t,2>){nidx1, nidx2}) = epsilons_[i];
      sigma_matrix_.at((std::array<std::size_t,2>){nidx1, nidx2}) = sigmas_[i];
      epsilon_matrix_.at((std::array<std::size_t,2>){nidx2, nidx1}) = epsilons_[i];
      sigma_matrix_.at((std::array<std::size_t,2>){nidx2, nidx1}) = sigmas_[i];
    }
    return;
  }
  virtual double compute_potential(int idx1, int idx2){
    Vec3<double> x1,x2;
    double sigma_effective, epsilon_effective;
    int midx1, midx2;
    auto mptr1 = index_map_.find(box->atoms[idx1].name);
    auto mptr2 = index_map_.find(box->atoms[idx2].name);
    FANCY_ASSERT(mptr1 != index_map_.end(), "Calc_SR_LJAtt() index not found in map");
    FANCY_ASSERT(mptr2 != index_map_.end(), "Calc_SR_LJAtt() index not found in map");
    midx1 = mptr1->second;
    midx2 = mptr2->second;
    sigma_effective = sigma_matrix_.at((std::array<int,2>){midx1, midx2});
    epsilon_effective = epsilon_matrix_.at((std::array<int,2>){midx1, midx2});
    x1=box->atoms[idx1].x;
    x2=box->atoms[idx2].x;
    getNearestImage3D(x2,x1,box_size_);
    double r = getDistanceNoPBC(x1,x2);
    if(r > cutoff_) return 0.0;
    if(r < sigma16_prefactor_*sigma_effective) return -epsilon_effective;
    return LJ_6_12(r, epsilon_effective, sigma_effective);
  }
protected:
  std::vector<double> sigmas_, epsilons_;
  std::vector<std::string> atom_names_, atom_names2_;
  std::map<std::string, std::size_t> index_map_;
  Matrix<double, 2> sigma_matrix_, epsilon_matrix_;
  double rcut_;
  bool shift_ = 1;
  double sigma16_prefactor_;
};