#pragma once
#include <unordered_map>
#include "Calc_1D_Density.hpp"
#include "Calc_HBonds.hpp"
#include "../probevolumes/ProbeVolume.hpp"
#include "../../tools/cellgrid.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/Matrix.hpp"
//scheme clusters, omitting grid densities below thresh
//each cluster is averaged, with the sum needing t o

class Calc_1D_HBond : public Calc_1D_Density{
  public:
    Calc_1D_HBond(InputPack& input):Calc_1D_Density{input}{
      FANCY_ASSERT(box->hasNamedAtoms, "Charge mapping scheme requires atom names to function");
      std::vector<std::string> atom_names;
      std::vector<double> charges, epsilons, sigmas;
      input.params().readVector("atom_names", ParameterPack::KeyType::Required,  atom_names);
      input.params().readVector("charges", ParameterPack::KeyType::Required, charges);
      FANCY_ASSERT(atom_names.size() == charges.size(), "Need one charge per atom name.");
      for(int i = 0; i < atom_names.size(); i++){
        charge_map_.insert(std::pair<std::string, double>(atom_names[i], charges[i]));
      }

      auto indices = atom_group_->getIndices();
      for(auto& index : indices){
        auto& atom = box->atoms[index];
        FANCY_ASSERT(charge_map_.find(atom.name) != charge_map_.end(), "Failed to find entry in map for one of the atoms in the group.");
      }
      //have ensured that every atom in the atom group has a corresponding charge in the map
      input.params().readNumber("cutoff", KeyType::Required, cutoff_);
      std::string calc_name_;
      input.params().readString("calc", KeyType::Required, calc_name_);
      Calculation* calc = input.findCalculation(calc_name_);
      FANCY_ASSERT(calc != 0, "Failed to find calculation in Calc_DensityField_HBond() : " + name_);
      calc_ = dynamic_cast<Calc_HBonds*>(calc);
      angle_density_.resize(grid_density_.size(), 0.0), average_angle_density_.resize(grid_density_.size(), 0.0);
      energy_density_.resize(grid_density_.size(), 0.0), average_energy_density_.resize(grid_density_.size(), 0.0);
      costheta_density_.resize(grid_density_.size(), 0.0), average_costheta_density_.resize(grid_density_.size(), 0.0);
      std::string pvname;
      input.params().readString("probevolume", KeyType::Required, pvname);
      pv_ = input.findProbeVolume(pvname);
      FANCY_ASSERT(pv_ != 0, "Failed to find probevolume in Calc_1D_Hbond() : " + name_ );
      return;
    }
    void putInBinWeighted(Vec3<double> pos, double weight, Vec<double>& grid_){
      int idx = getBin(pos[ dim_ ]);
      grid_[idx] += weight;
      return;
    }
    virtual void update(){
      Calc_1D_Density::update();
      for(int i = 0; i < grid_density_.size(); i++){
        energy_density_[i] = 0.0;
        angle_density_[i] = 0.0;
        costheta_density_[i] = 0.0;
      }     
    }
  void add_gaussian_weighted(double x_in, double weight, Vec<double>& grid_)
  {
      double x = x_in;
      int lxmin = floor((x-2*sigma_)/grid_spacing_);
      int lxmax = ceil((x+2*sigma_)/grid_spacing_); 
      #pragma omp parallel for
      for(int ix = lxmin; ix <= lxmax; ix++)
      {
        int idx;
        idx = ix;
        if(idx >= grid_.size()) idx -= grid_.size();
        else if(idx < 0) idx += grid_.size();
        double xmin, xmax;
        xmin = ix * grid_spacing_;
        xmax = xmin + grid_spacing_;
        grid_[idx] += h_x(x, xmin, xmax, sigma_, 2.0*sigma_)*weight;
      }
      return;
  }
    virtual void calculate(){
      if(!doCalculate()) return;
      if(!calc_->hasCalculated()) calc_->calculate();
      auto hbonds = calc_->getHBonds();
      
      #pragma omp parallel for
      for(const auto& hbond : hbonds){
        auto position = hbond.x;
        auto inProbe = pv_->compute(position);
        std::string name1 = box->atoms[hbond.indices[0]].name, name2 = box->atoms[hbond.indices[1]].name;
        double q1 = charge_map_.find(name1)->second;
        double q2 = charge_map_.find(name2)->second;
        if(coarseGrain){
          add_gaussian_weighted(position[dim_], inProbe, grid_density_);
          add_gaussian_weighted(position[dim_], inProbe*Usr_test(hbond.r, q1*q2), energy_density_);
          add_gaussian_weighted(position[dim_], inProbe*hbond.angle, angle_density_);
          add_gaussian_weighted(position[dim_], inProbe * hbond.v[dim_] * 1.0/hbond.r, costheta_density_);
        }
        else{
          putInBinWeighted(position, inProbe, grid_density_);
          putInBinWeighted(position, inProbe*Usr_test(hbond.r, q1*q2), energy_density_);
          putInBinWeighted(position, inProbe*hbond.angle, angle_density_);
          putInBinWeighted(position, inProbe * hbond.v[dim_] * 1.0/hbond.r, costheta_density_);
        }
      }
      for(int i = 0; i < grid_density_.size(); i++){
        average_grid_density_[i] += grid_density_[i];
        average_angle_density_[i] += angle_density_[i];
        average_energy_density_[i] += energy_density_[i];
        average_costheta_density_[i] += costheta_density_[i];
      }
      average_grid_spacing_ += grid_spacing_;
      frame_counter_++;
    }
  void finalOutput(){
    average_grid_density_ = average_grid_density_* (1.0/(double)frame_counter_);
    average_energy_density_ = average_energy_density_* (1.0/(double)frame_counter_);
    average_angle_density_ = average_angle_density_* (1.0/(double)frame_counter_);
    average_costheta_density_ = average_costheta_density_* (1.0/(double)frame_counter_);
    for(int i = 0; i < average_grid_density_.size(); i++){
      average_energy_density_[i] /= average_grid_density_[i];
      average_angle_density_[i] /= average_grid_density_[i];
      average_costheta_density_[i] /= average_grid_density_[i];
    }
    average_grid_spacing_*= 1.0/frame_counter_;
    std::ofstream ofile(base_ + "_1D.txt");
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for 1D density calculation.");
    for(int i = 0; i < average_grid_density_.size(); i++){
      ofile << (i+0.5)*average_grid_spacing_ << "     " << average_grid_density_[i] << "   "
            << average_angle_density_[i]  << "   "
            << average_energy_density_[i]  << "   "
            << average_costheta_density_[i] << "   "
            << std::endl;
    }
    ofile.close();

    return;
  }


  protected:
    //usr between nearby atom of charge q and test charge
    inline double Usr_test(double  r,  double q){
      if(r < 0.1) r = 0.1;
      return (q/r)*(138.935458*0.5);
    }
    std::map<std::string, double> charge_map_;
    std::map<std::string, double> eps_map_, sigma_map_;
    double cutoff_;
    Calc_HBonds* calc_;
    ProbeVolume* pv_;
    Vec<double> angle_density_, average_angle_density_;
    Vec<double> energy_density_, average_energy_density_;
    Vec<double> costheta_density_, average_costheta_density_;
};