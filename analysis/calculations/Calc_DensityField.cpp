#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/smearfunctions.hpp"
Calc_DensityField::Calc_DensityField(InputPack& input):Calculation{input}{
  std::string agname;
  span_.fill(0);
  sigma_ = 0.0;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find atom group.");
  coarseGrain_ = 0;
  input.params().readFlag("smear", KeyType::Optional, coarseGrain_);
  if(coarseGrain_){
    input.params().readNumber("sigma", KeyType::Required, sigma_);
    FANCY_ASSERT(sigma_ > 0.0, "sigma < 0 is invalid");
  }
  std::vector<int> npoints;
  std::vector<double> box;
  input.params().readVector("npoints", KeyType::Required, npoints);
  for(int i = 0; i < 3; i++){
    npoints_[i] = npoints[i];
  }
  FANCY_ASSERT(npoints.size() == 3, "Improper size for npoints. Needs 3: nx, ny, nz");
  hasBoxVec_ = input.params().readVector("box_vector", KeyType::Optional, box);
  if(hasBoxVec_){
    FANCY_ASSERT(box.size() == 6, "Improper size for box. Needs 6: min_x, min_y, min_z, max_x, max_y, max_z");
    for(int i = 0; i<3; i++){
      minx_[i] = box[i];
      maxx_[i] = box[i+3];
      FANCY_ASSERT(minx_[i] < maxx_[i], "box has invalid dimensions, min_x < max_x");
      gridspacing_[i] = (maxx_[i] - minx_[i]) / (double)npoints_[i];
      if(coarseGrain_){
        span_[i] = ceil(2.0*sigma_/gridspacing_[i]);
      }
    }
  }
  gridvals_.resize(npoints_[0]*npoints_[1]*npoints_[2], 0.0);
  avggridvals_.resize(npoints_[0]*npoints_[1]*npoints_[2], 0.0);
  avggridspacing_.fill(0.0);
  nframes_ = 0;
  return;
}
Vec3<int> Calc_DensityField::getIndex(const Vec3<double>& pos){
  Vec3<int> idx3;
  for(int i = 0; i < 3; i++){
    idx3[i] = floor((pos[i] - minx_[i]) / gridspacing_[i]), npoints_[i];
  }
  return idx3;
}
bool Calc_DensityField::getIndexNoWrap1D(double& pos, int dim, int& ret){
  Vec3<int> idx_ref;
  auto newpos = pos;
  wrapNumber(newpos, box_size_[dim]);
  ret = floor((newpos - minx_[dim]) / gridspacing_[dim]);
  if(idx_ref[dim] < 0 - span_[dim] || idx_ref[dim] > npoints_[dim] + span_[dim]) return 0;
  return 1;
}

double Calc_DensityField::Calc_DensityField::getGaussian(const Vec3<double>& x, const Vec3<int>& index){
  Vec3<double> xmin, xmax;
  double eval = 1.0;
  for(int i = 0; i < 3; i++){
    double x_pbc = getNearestImage1D(x[i], index[i]*gridspacing_[i] + minx_[i], box_size_[i]);
    xmin[i] = minx_[i] + index[i]*gridspacing_[i];
    xmax[i] = xmin[i] + gridspacing_[i];
    eval *= h_x(x_pbc, xmin[i], xmax[i], sigma_, 2.0*sigma_);
  }
  return eval;
}
void Calc_DensityField::Calc_DensityField::calculate(){
  if(!doCalculate()) return;
  for(auto& value : gridvals_){
    value = 0.0;
  }
  auto indices = atom_group_->getIndices();
  for(auto index : indices){
    auto position = box->atoms[index].x;
    placeInsideBox(position, box_size_);
    bool flag = 0;
    for(int i = 0; i < 3; i++){
      if(position[i] < minx_[i] - 2.0*sigma_ || position[i] > maxx_[i] + 2.0*sigma_){
        flag = 1;
        break;
      }
      }
    if(flag) continue;
    auto idx_ref = getIndex(position);
    //will not use periodic boundaries if using a subvolume
    if(coarseGrain_){
      Vec3<int> idx_temp;
      if(hasBoxVec_){
        for(int i = -span_[0]; i <= span_[0]; i++){
          idx_temp[0] = idx_ref[0] + i;
          if(idx_temp[0] < 0 || idx_temp[0] >= npoints_[0]) continue;
          for(int j = -span_[1]; j <= span_[1]; j++){
            idx_temp[1] = idx_ref[1] + j;
            if(idx_temp[1] < 0 || idx_temp[1] >= npoints_[1]) continue;
            for(int k = -span_[2]; k <= span_[2]; k++){
              idx_temp[2] = idx_ref[2] + k;
              if(idx_temp[2] < 0 || idx_temp[2] >= npoints_[2]) continue;
              gridvals_[_map31(idx_temp)] += getGaussian(position, idx_temp);
            }
          }
        }
      }
      else{
        for(int i = -span_[0]; i <= span_[0]; i++){
          idx_temp[0] = wrapIndex(idx_ref[0] + i, npoints_[0]);
          if(idx_temp[0] < 0 || idx_temp[0] >= npoints_[0]) continue;
          for(int j = -span_[1]; j <= span_[1]; j++){
            idx_temp[1] = wrapIndex(idx_ref[1] + j, npoints_[1]);
            for(int k = -span_[2]; k <= span_[2]; k++){
              idx_temp[2] = wrapIndex(idx_ref[2], npoints_[2]);
              gridvals_[_map31(idx_temp)] += getGaussian(position, idx_temp);
            }
          }
        }        
      }
    }
    else{
      gridvals_[_map31(idx_ref)] += 1.0;
    }
  }
  avggridvals_ = avggridvals_ + gridvals_;
  avggridspacing_ = avggridspacing_ + gridspacing_;
  nframes_++;
}
void Calc_DensityField::update(){
  Calculation::update();
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  if(!hasBoxVec_){
    for(int i = 0; i < 3; i++){
      minx_[i] = 0;
      maxx_[i] = box->boxvec[i][i];
      gridspacing_[i] = (maxx_[i] - minx_[i]) / (double)npoints_[i];
      if(coarseGrain_){
        span_[i] = ceil(2.0*sigma_/gridspacing_[i]);
      }
    }
  }
  return;
}
std::string Calc_DensityField::printConsoleReport(){
  return "";
}
void Calc_DensityField::finalOutput(){
  avggridvals_ = (1.0/(double)nframes_) * avggridvals_;
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
  return;
}