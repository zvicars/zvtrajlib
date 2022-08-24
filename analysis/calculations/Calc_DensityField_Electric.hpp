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
      input.params().readNumber("beta", KeyType::Required, beta_);
      //alpha is set based on the peak width being 0.5/rc
      sigma_ = cutoff_/2.0;
      return;
    }
    virtual void calculate(){
      if(!doCalculate()) return;
      for(auto& value : gridvals_){
        value = 0.0;
      }
      for(int i = 0; i < 3; i++){
        span_[i] = ceil(cutoff_ / gridspacing_[i]);
      }
      auto atomIndices = atom_group_->getIndices();
      #pragma omp parallel for
      for(auto idx : atomIndices){
        auto atom = box->atoms[idx];
        auto pos = atom.x;
        Vec3<double> min_pos;
        Vec3<double> max_pos; 
        for(int i = 0; i < 3; i++){
          min_pos[i] = std::fmax((pos[i] - minx_[i]) / gridspacing_[i] - span_[i], 0);
          max_pos[i] = std::fmin((pos[i] - minx_[i]) / gridspacing_[i] + span_[i], npoints_[i]-1);
        }
        Vec3<double> r1;
        for(int i = min_pos[0]; i <= max_pos[0]; i++){
          for(int j = min_pos[1]; j <= max_pos[1]; j++){
            for(int k = min_pos[2]; k <= max_pos[2]; k++){
              r1[0] = i*gridspacing_[0] + minx_[0];
              r1[1] = j*gridspacing_[1] + minx_[1];
              r1[2] = k*gridspacing_[2] + minx_[2];
              Vec3<int> idx3d = {i,j,k};
              auto r = getDistance(r1, pos, box_size_);
              if(r > cutoff_) continue;
              auto idx1d = _map31(idx3d);
              gridvals_[idx1d] += Usr_test(r, charge_map_.find(atom.name)->second);
            }
          }
        }
      }
      gridvals_ = gridvals_ * (138.935458*0.5);
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
    std::map<std::string, double> charge_map_;
    double cutoff_, sigma_, alpha_, beta_;
};