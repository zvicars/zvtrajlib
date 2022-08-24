#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
#include "../../interface/interface.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/cellgrid.hpp"
#include <unordered_map>
#include <cstdlib>
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
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
      for(int i = 0; i < 3; i++){
        span_[i] = ceil(cutoff_ / gridspacing_[i]);
      }
      auto atomIndices = atom_group_->getIndices();
      //fast but doesn't work across pbc's
      //pragma omp parallel for
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
              gridvals_[idx1d] += LJ_single(r, sigma_map_.find(atom.name)->second, eps_map_.find(atom.name)->second);
            }
          }
        }
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
    std::map<std::string, double> eps_map_, sigma_map_;
    double cutoff_, sigma_, alpha_, beta_;
    double max_spacing_;
};