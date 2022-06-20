#include "Calc_DensityField.hpp"
#include "Eigen/Eigen"
#include "../../tools/stlmath.hpp"
#include "../../interface/interface.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../helper/make_histogram.hpp"
#include <unordered_map>
//special angle calculation to complete the tetrahedron for water

class Calc_DensityFieldWater : public Calc_DensityField{
  public:
    Calc_DensityFieldWater(InputPack& input):Calc_DensityField{input}{
      FANCY_ASSERT(coarseGrain_ == 0, "Coarse graining is not supported for this derived class");
      std::string ag2, ag3;
      input.params().readString("atom_group2", ParameterPack::KeyType::Required, ag2);
      atom_group2_ = input.findAtomGroup(ag2);
      FANCY_ASSERT(atom_group2_ != 0, "Failed to find second atom group in Calc_DensityFieldWater()");
      input.params().readString("atom_group3", ParameterPack::KeyType::Required, ag3);
      atom_group3_ = input.findAtomGroup(ag3);
      FANCY_ASSERT(atom_group3_ != 0, "Failed to find third atom group in Calc_DensityFieldWater()");
      std::vector<double> normal_vector;
      input.params().readVector("normal_vector", ParameterPack::KeyType::Required, normal_vector);
      FANCY_ASSERT(normal_vector.size() == 3, "Improper size for normal vector in Calc_DensityFieldAngle()");
      normal_vector_ = vec2Array<double,3>(normal_vector);
      normal_vector_ = normal_vector_ * (1.0/norm2(normal_vector_));
      gridsum1_.resize(npoints_[0]*npoints_[1]*npoints_[2], 0.0);
      gridsum2_.resize(npoints_[0]*npoints_[1]*npoints_[2], 0.0);
      gridsum3_.resize(npoints_[0]*npoints_[1]*npoints_[2], 0.0);

      input.params().readFlag("histogram", KeyType::Optional, doHistogram);;
      bool found_min = input.params().readNumber("min_bin", KeyType::Optional, min_bin_);
      if(found_min) forceMin = 1;
      bool found_max = input.params().readNumber("max_bin", KeyType::Optional, max_bin_);
      if(found_max) forceMax = 1;
      bool found_bs = input.params().readNumber("bin_size", KeyType::Optional, bin_size_);
      if(found_bs) forceBS = 1;
      if(found_min || found_max || found_bs) doHistogram = 1;

      return;
    }
    virtual void calculate(){
      Calc_DensityField::calculate();
      if(!doCalculateChild()) return;
      FANCY_ASSERT(atom_group_->getIndexCount() == atom_group2_->getIndexCount() && atom_group_->getIndexCount() == atom_group3_->getIndexCount(), "Atom groups are not the same size in Calc_DensityFieldAngle().");
      auto& idx1 = atom_group_->getIndices();
      auto& idx2 = atom_group2_->getIndices();
      auto& idx3 = atom_group3_->getIndices();
      for(int i = 0; i < idx1.size(); i++){
        int index1 = idx1[i];
        int index2 = idx2[i];
        int index3 = idx3[i];
        //vectors are OW-bisecting, OW-HW1, OW-HW2, OW-tet1, OW-tet2
        std::array<double,3> v1,v2,v3,v4,v5;
        auto x_ow = box->atoms[index1].x;
        if(!isInBox(x_ow)) continue; 
        placeInsideBox(x_ow, box_size_);
        auto x_hw1 = box->atoms[index2].x;
        auto x_hw2 = box->atoms[index3].x;
        getNearestImage3D(x_hw1, x_ow, box_size_);
        getNearestImage3D(x_hw2, x_ow, box_size_);
        buildVectors(x_ow, x_hw1, x_hw2, v1, v2, v3, v4, v5 );
        double angle1 = dot(v1, normal_vector_);
        double angle2 = dot(v2, normal_vector_);
        double angle3 = dot(v3, normal_vector_);
        double angle4 = dot(v4, normal_vector_);
        double angle5 = dot(v5, normal_vector_);
        auto idx_ref = getIndex(x_ow);
        gridsum1_[_map31(idx_ref)] += angle1; 
        gridsum2_[_map31(idx_ref)] += 0.5*(angle2 + angle3);
        gridsum3_[_map31(idx_ref)] += 0.5*(angle4 + angle5);
        ts1_.push_back(angle1);
        ts2_.push_back(angle2);
        ts2_.push_back(angle3);
        ts3_.push_back(angle4);
        ts3_.push_back(angle5); 
      }
      return;
    }

    virtual void finalOutput(){
      for(int i = 0; i < avggridvals_.size(); i++){
        if(avggridvals_[i] != 0.0){
          gridsum1_[i]    *= 1.0/avggridvals_[i];
          gridsum2_[i]    *= 1.0/avggridvals_[i];
          gridsum3_[i]    *= 1.0/avggridvals_[i];
        }
        else{
          gridsum1_[i]    = 0.0;
          gridsum2_[i]    = 0.0;
          gridsum3_[i]    = 0.0;
        }
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
            ofile << i << "   " << j << "   " << k << "   " << avggridvals_[idx] << "   " 
                  << gridsum1_[idx] << "   " << gridsum2_[idx] << "   " << gridsum3_[idx] << "\n";
          }
        }
      }
      ofile.close();

      double mean1 = mean(ts1_);
      double var1 = var(ts1_, mean1);
      double mean2 = mean(ts2_);
      double var2 = var(ts2_, mean2);
      double mean3 = mean(ts3_);
      double var3 = var(ts3_, mean3);


      ofile.open(base_ + "_statistics.txt");
      FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
      ofile << "#Output file for DensityField_Water calculation with name \'" << name_ << "\'\n";
      ofile << "#Atom group: " << atom_group_->getName() << "\n";
      ofile << "#Average: " << mean1 << "   " << mean2  << "   " << mean3 << "\n";
      ofile << "#Variance: " << var1 << "   " << var2  << "   " << var3 << "\n";
      if(doHistogram){
          std::vector<double> x_vals;
          std::vector<int> y_vals;
          ofile << "#Histogram: cos1     count\n";
          makeHistogram(ts1_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
          for(std::size_t i = 0; i < x_vals.size(); i++){
            ofile << x_vals[i] << "   " << y_vals[i] << "\n";
          }
          ofile << "#Histogram: cos2     count\n";
          makeHistogram(ts2_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
          for(std::size_t i = 0; i < x_vals.size(); i++){
            ofile << x_vals[i] << "   " << y_vals[i] << "\n";
          }
          ofile << "#Histogram: cos3     count\n";
          makeHistogram(ts3_, min_bin_, max_bin_, bin_size_, forceMin, forceMax, forceBS, x_vals, y_vals);
          for(std::size_t i = 0; i < x_vals.size(); i++){
            ofile << x_vals[i] << "   " << y_vals[i] << "\n";
          }    
      }
      ofile.close();
      return;
    }
  protected:
    void buildVectors(const Vec3<double>& x_ow, const Vec3<double>& x_hw1, const Vec3<double>& x_hw2, 
    Vec3<double>& v1, Vec3<double>& v2, Vec3<double>& v3, Vec3<double>& v4, Vec3<double>& v5){
        //vector 1 goes from the OW atom to the midpoint of the HW1 and HW2 atoms
        Vec3<double> x_mid;
        x_mid = (x_hw1 + x_hw2) * 0.5;
        v1 = (1.0/norm2(x_ow - x_mid)) * (x_ow - x_mid);
        v2 = (1.0/norm2(x_hw1 - x_ow)) * (x_hw1 - x_ow);
        v3 = (1.0/norm2(x_hw1 - x_ow)) * (x_hw2 - x_ow);
        //hard part is building tetrahedron...
        //bisector, v1, should bisect tetrahedron, v4, v5 should lie on plane 90 deg offset from v2,v3
        Eigen::Vector3d ev1, ev2, ev3, ev4, ev5;
        ev1 << v1[0], v1[1], v1[2];
        ev2 << v2[0], v2[1], v2[2];
        ev3 << v3[0], v3[1], v3[2];
        //plane 1: ow hw1 hw2
        auto en1 = ev2.cross(ev3);
        en1.normalize();
        //plane 2 can be found by rotating 90 degrees about v1
        auto m = Eigen::AngleAxisd(M_PI_2, ev1);
        auto en2 = m*en1;
        //vectors 4,5 can be found by the +- 109.5/2 rotation about the en2 normal
        #define ROT_ANGLE 0.955567765
        m = Eigen::AngleAxisd(ROT_ANGLE, en2);
        ev4 = m*ev1;
        m = Eigen::AngleAxisd(-ROT_ANGLE, en2);
        ev5 = m*ev1;
        for(int i = 0; i < 3; i++){
          v4[i] = ev4[i];
          v5[i] = ev5[i];
        }
        v1 = v1 / norm2(v1);
        v2 = v2 / norm2(v2);
        v3 = v3 / norm2(v3);
        v4 = v4 / norm2(v4);
        v5 = v5 / norm2(v5);
      return;
    }
    AtomGroup* atom_group2_, * atom_group3_;
    std::vector<double> gridsum1_, gridsum2_, gridsum3_;
    std::vector<double> ts1_, ts2_, ts3_;
    std::array<double,3> normal_vector_;
    double bin_size_, min_bin_, max_bin_; //standard histogram output options, min/max bin can be found dynamically too
    bool doHistogram=0, forceMin=0, forceMax=0, forceBS=0;
};