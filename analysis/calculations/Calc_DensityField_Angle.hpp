#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
#include "../../interface/interface.hpp"
#include "../../tools/pbcfunctions.hpp"
#include <unordered_map>
//scheme clusters, omitting grid densities below thresh
//each cluster is averaged, with the sum needing t o

class Calc_DensityFieldAngle : public Calc_DensityField{
  public:
    Calc_DensityFieldAngle(InputPack& input):Calc_DensityField{input}{
      FANCY_ASSERT(box->hasNamedAtoms, "Charge mapping scheme requires atom names to function");
      std::string ag2_name;
      input.params().readString("atom_group2", ParameterPack::KeyType::Required,  ag2_name);
      atom_group2_ = input.findAtomGroup(ag2_name);
      FANCY_ASSERT(atom_group2_ != 0, "Failed to find second atom group in Calc_DensityFieldAngle()");
      std::vector<double> normal_vector;
      input.params().readVector("normal_vector", ParameterPack::KeyType::Required, normal_vector);
      FANCY_ASSERT(normal_vector.size() == 3, "Improper size for normal vector in Calc_DensityFieldAngle()");
      normal_vector_ = vec2Array<double,3>(normal_vector);
      normal_vector_ = normal_vector_ * (1.0/norm2(normal_vector_));
      gridcounts_.resize(npoints_[0]*npoints_[1]*npoints_[2], 0.0);
      return;
    }
    virtual void calculate(){
      if(!doCalculate()) return;
      for(auto& value : gridvals_){
        value = 0.0;
      }
      FANCY_ASSERT(atom_group_->getIndexCount() == atom_group2_->getIndexCount(), "Atom groups are not the same size in Calc_DensityFieldAngle().");
      auto& idx1 = atom_group_->getIndices();
      auto& idx2 = atom_group2_->getIndices();
      for(int i = 0; i < idx1.size(); i++){
        int index1 = idx1[i];
        int index2 = idx2[i];
        auto& x1 = box->atoms[index1].x;
        auto x2 = box->atoms[index2].x;
        getNearestImage3D(x2, x1, box_size_);
        auto x21 = (x2-x1) * (1.0/norm2(x2-x1));
        double angle = dot(x21, normal_vector_);
        auto idx_ref = getIndex(x1);
        gridvals_[_map31(idx_ref)] += angle;
        gridcounts_[_map31(idx_ref)] += 1; 
      }
      nframes_++;
      avggridvals_ = avggridvals_ + gridvals_;
      avggridspacing_ = avggridspacing_ + gridspacing_;
      return;
    }

    virtual void finalOutput(){
      for(int i = 0; i < avggridvals_.size(); i++){
        avggridvals_[i] *= 1.0/gridcounts_[i];
      }
      avggridspacing_ = (1.0/(double)nframes_)  * avggridspacing_;
      std::ofstream ofile(base_ + "_DensityField.txt");
      FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
      for(int i = 0; i < npoints_[0]; i++){
        for(int j = 0; j < npoints_[1]; j++){
          for(int k = 0; k < npoints_[2]; k++)
          {
            Vec3<int> pos = {i,j,k};
            int idx = _map31(pos);
            ofile << i << "   " << j << "   " << k << "   " << avggridvals_[idx] << "\n";
          }
        }
      }
      ofile.close();
      return;
    }
  protected:
    AtomGroup* atom_group2_;
    std::vector<double> gridcounts_;
    std::array<double,3> normal_vector_;
};