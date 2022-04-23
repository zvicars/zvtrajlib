#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
#include "../../interface/interface.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/cellgrid.hpp"
#include <unordered_map>
//scheme clusters, omitting grid densities below thresh
//each cluster is averaged, with the sum needing t o

class Calc_DensityFieldElectric : public Calc_DensityField{
  public:
    Calc_DensityFieldElectric(InputPack& input):Calc_DensityField{input}{
      FANCY_ASSERT(box->hasNamedAtoms, "Charge mapping scheme requires atom names to function");
      std::vector<std::string> atom_names;
      std::vector<double> charges;
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
      //alpha is set based on the peak width being 0.5/rc
      sigma_ = cutoff_/2.0;
      alpha_ = 2.0/std::pow((2*sigma_*sigma_),2);
      beta_ = std::sqrt(alpha_);
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
        gridvals_[i] = Usr_test_sum(realpos_);
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
    inline double Usr_test(double  r,  double q){
      if(r < 0.1*norm2(gridspacing_)) r = 0.1*norm2(gridspacing_);
      return (q/r)*erfc(beta_*r);
    }
    inline void resetCellGrid(){
      cell_grid_.reset(cutoff_, box_size_);
      auto& indices = atom_group_->getIndices();
      for(auto index : indices){
        cell_grid_.addIndexToGrid(index, box->atoms[index].x);
      }
    }
    //usr summed over all atoms for a given position with a test charge of 1
    inline double Usr_test_sum(std::array<double,3> r1){
      double sum = 0;
      auto nearbyAtomIndices = cell_grid_.getNearbyIndices(r1);
      for(auto idx : nearbyAtomIndices){
        auto& atom = box->atoms[idx];
        auto r = getDistance(r1, atom.x, box_size_);
        sum += Usr_test(r, charge_map_.find(atom.name)->second);
      }
      return sum; 
    }
    std::map<std::string, double> charge_map_;
    double cutoff_, sigma_, alpha_, beta_;
    CellGrid cell_grid_;
};