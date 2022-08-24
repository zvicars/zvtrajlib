#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
#include "../../interface/interface.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/cellgrid.hpp"
#include <unordered_map>

class Calc_DensityFieldH2OMin : public Calc_DensityField{
  public:
    Calc_DensityFieldH2OMin(InputPack& input):Calc_DensityField{input}{
      FANCY_ASSERT(box->hasNamedAtoms, "Charge mapping scheme requires atom names to function");
      std::vector<std::string> atom_names;
      std::vector<double> charges, epsilons, sigmas;
      input.params().readVector("atom_names", ParameterPack::KeyType::Required,  atom_names);
      input.params().readVector("charges", ParameterPack::KeyType::Required, charges);
      input.params().readVector("epsilons", ParameterPack::KeyType::Required, epsilons);
      input.params().readVector("sigmas", ParameterPack::KeyType::Required, sigmas);

      FANCY_ASSERT(atom_names.size() == charges.size(), "Need one charge per atom name.");
      for(int i = 0; i < atom_names.size(); i++){
        charge_map_.insert(std::pair<std::string, double>(atom_names[i], charges[i]));
      }

      FANCY_ASSERT(atom_names.size() == sigmas.size(), "Need one sigma per atom name.");
      FANCY_ASSERT(atom_names.size() == epsilons.size(), "Need one epsilon per atom name.");
      for(int i = 0; i < atom_names.size(); i++){
        eps_map_.insert(std::pair<std::string, double>(atom_names[i], epsilons[i]));
        sigma_map_.insert(std::pair<std::string, double>(atom_names[i], sigmas[i]));
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

      normals_.resize(npoints_[0]*npoints_[1]*npoints_[2], {0,0,0});
      avgNormals_.resize(npoints_[0]*npoints_[1]*npoints_[2], {0,0,0});
      return;
    }
    virtual void calculate();
    virtual void finalOutput();
  protected:
    virtual double minimize_ff(Vec3<double> site_position, const std::vector<int>& indices, Vec3<double>& eulerStep);
    //usr between nearby atom of charge q and test charge
    inline double Usr_test(double  r,  double q){
      if(r < 0.1*norm2(gridspacing_)) r = 0.1*norm2(gridspacing_);
      return (q/r)*erfc(beta_*r);
    }
    std::vector<Vec3<double>> normals_;
    std::vector<Vec3<double>> avgNormals_;
    Vec3<double> current_euler_angles_;
    std::map<std::string, double> charge_map_, eps_map_, sigma_map_;
    double cutoff_, sigma_, alpha_, beta_;
};