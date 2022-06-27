#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
#include "../../interface/interface.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/cellgrid.hpp"
#include <unordered_map>
//scheme clusters, omitting grid densities below thresh
//each cluster is averaged, with the sum needing t o

class Calc_DensityFieldLJ : public Calc_DensityField{
  public:
    Calc_DensityFieldLJ(InputPack& input):Calc_DensityField{input}{
      FANCY_ASSERT(box->hasNamedAtoms, "Charge mapping scheme requires atom names to function");
      std::vector<std::string> atom_names;
      std::vector<double> epsilons, sigmas;
      input.params().readVector("atom_names", ParameterPack::KeyType::Required,  atom_names);
      input.params().readVector("epsilons", ParameterPack::KeyType::Required, epsilons);
      input.params().readVector("sigmas", ParameterPack::KeyType::Required, sigmas);
      FANCY_ASSERT(atom_names.size() == sigmas.size(), "Need one sigma per atom name.");
      FANCY_ASSERT(atom_names.size() == epsilons.size(), "Need one epsilon per atom name.");
      for(int i = 0; i < atom_names.size(); i++){
        eps_map_.insert(std::pair<std::string, double>(atom_names[i], epsilons[i]));
        sigma_map_.insert(std::pair<std::string, double>(atom_names[i], sigmas[i]));
      }

      auto indices = atom_group_->getIndices();
      for(auto& index : indices){
        auto& atom = box->atoms[index];
        FANCY_ASSERT(eps_map_.find(atom.name) != eps_map_.end(), "Failed to find entry in map for one of the atoms in the group.");
        FANCY_ASSERT(sigma_map_.find(atom.name) != sigma_map_.end(), "Failed to find entry in map for one of the atoms in the group.");
      }
      //have ensured that every atom in the atom group has a corresponding charge in the map
      input.params().readNumber("cutoff", KeyType::Required, cutoff_);
      return;
    }
    virtual void calculate(){
      if(!doCalculate()) return;
      for(auto& value : gridvals_){
        value = 0.0;
      }
      resetCellGrid();
      #pragma omp parallel for
      for(int i = 0; i < gridvals_.size(); i++){
        auto idx3d = _map13(i);
        std::array<double, 3> realpos_;
        for(int j = 0; j < 3; j++){
          realpos_[j] = idx3d[j]*gridspacing_[j] + minx_[j];
        }
        gridvals_[i] = LJ_sum(realpos_);
      }
      nframes_++;
      avggridvals_ = avggridvals_ + gridvals_;
      avggridspacing_ = avggridspacing_ + gridspacing_;
      return;
    }
    virtual void finalOutput(){
      Calc_DensityField::finalOutput();
      return;
    }
  protected:
    //usr between nearby atom of charge q and test charge
    inline double LJ_single(double  r,  double sigma, double epsilon){
      if(r < 0.1*norm2(gridspacing_)) r = 0.1*norm2(gridspacing_);
      double r6 = std::pow(sigma/r,6.0);
      double r12 = r6*r6;
      return -4.0*epsilon*(r6 - r12);
    }
    inline void resetCellGrid(){
      cell_grid_.reset(cutoff_, box_size_);
      auto& indices = atom_group_->getIndices();
      for(auto index : indices){
        cell_grid_.addIndexToGrid(index, box->atoms[index].x);
      }
    }
    //usr summed over all atoms for a given position with a test charge of 1
    inline double LJ_sum(std::array<double,3> r1){
      double sum = 0;
      auto nearbyAtomIndices = cell_grid_.getNearbyIndices(r1);
      for(auto idx : nearbyAtomIndices){
        auto& atom = box->atoms[idx];
        auto r = getDistance(r1, atom.x, box_size_);
        sum += LJ_single(r, sigma_map_.find(atom.name)->second, eps_map_.find(atom.name)->second);
      }
      return sum; 
    }
    std::map<std::string, double> eps_map_, sigma_map_;
    double cutoff_, sigma_, alpha_, beta_;
    CellGrid cell_grid_;
};