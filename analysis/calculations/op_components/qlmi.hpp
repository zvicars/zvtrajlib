//calculates qbar values for a set of atoms
//needs location of atom and its neighbors
#pragma once
#include "spherical_harmonics.hpp"
#include "../../../tools/pbcfunctions.hpp"
#include "htildei.hpp"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class qlmi{
  public:
  qlmi(std::size_t order, double rmax, double sigma, double cutoff){
    l_ = order;
    v1 = new sphericalINDUSVolume({0,0,0}, rmax, sigma, cutoff);
    return;
  }
  ~qlmi(){
    delete v1;
  }
  inline int to_m(int i) const{
    return i - l_;
  }
  inline int to_idx(int m) const{
    return m + l_;
  }
  double qlmi_test_probe(Vec3<double> pos, Vec3<double> ref_pos){
    v1->updatePosition(ref_pos);
    v1->update_box_size({100,100,100});
    return v1->calculate(pos);
  }
  //returns number of nearest neighbors, updates ylm_sum_vec
  double compute(const Vec3<double>& atom_position, const Vec<Vec3<double> >& neighbor_positions, 
              Vec<std::complex<double> >& ylm_sum_vec, Vec<double>& htilde_neighbors) const{
    double eval=0.0;
    int nn = neighbor_positions.size();
    v1->updatePosition(atom_position);
    //if you want index for m, need to use i = m + l, i - l = m (i = l for m = 0)
    std::vector<std::vector<std::complex<double> > > ylm_vec(nn, std::vector<std::complex<double> >( 2*l_ + 1 ) );
    std::vector<double> htilde_vec(nn);
    for(int i = 0; i < nn; i++){
      Vec3<double> xn = neighbor_positions[i];
      getNearestImage3D(xn, atom_position, box_size_);
      htilde_vec[i] = v1->calculate(xn);
      xn = xn - atom_position; 
      double r = 0.0;
      for(int j = 0; j < 3; j++) r += xn[j]*xn[j];
      r = sqrt(r);
      double costheta = xn[2]/r;
      double phi = sgn(xn[1])*acos(xn[0]/sqrt(xn[0]*xn[0] + xn[1]*xn[1]));
      for(int m = 0; m <= l_; m++){
        ylm_vec[i][to_idx(m)] = Ylm(l_, m, costheta, phi);
        if(m != 0) ylm_vec[i][to_idx(-m)] = Ylm_invert( ylm_vec[i][to_idx(m)], m );
      }
    }
    double nni = 0.0;
    ylm_sum_vec.clear();
    ylm_sum_vec.resize(2*l_+1, std::complex<double>(0,0));
    for(int i = 0; i < nn; i++){
      nni += htilde_vec[i];
      for(int j = 0; j < 2*l_+1; j++ ){
        ylm_sum_vec[j] += ylm_vec[i][j]*htilde_vec[i];
      }
    }
    for(int i = 0; i < 2*l_+1; i++ ){
      if(nni > 0) ylm_sum_vec[i] *= 1.0/nni;
    }
    htilde_neighbors = htilde_vec;
    return nni;
  }
  void updateBoxSize(Vec3<double> box_size){
    v1->update_box_size(box_size);
    box_size_ = box_size;
  }
  protected:
  sphericalINDUSVolume* v1;
  Vec3<double> box_size_;
  bool prefactors_computed_ = 0;
  std::vector<double> prefactors_; //-l -l+1 .. 0 .. l-1 1, 2*l + 1 entries
  int l_;
};