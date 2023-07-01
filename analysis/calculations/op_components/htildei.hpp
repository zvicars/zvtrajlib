#pragma once
#include "../../../tools/Assert.hpp"
#include "../../../tools/stlmath.hpp"
#include "../../../interface/datatypes.hpp"
#include "../../../tools/smearfunctions.hpp"
#include "../../../tools/pbcfunctions.hpp"
class INDUSVolume{
public:
  INDUSVolume(){ };
  virtual double calculate(const Vec3<double>& x) const = 0;
  virtual void update_box_size(const Vec3<double>& box_size){
    box_size_ = box_size;
    return;
  }
protected:
  Vec3<double> box_size_;
};

class sphericalINDUSVolume : public INDUSVolume{
public:
  sphericalINDUSVolume(Vec3<double> position, double rmax, double sigma, double cutoff){
    position_ = position;
    sigma_ = sigma;
    rmax_ = rmax;
    rc_ = cutoff;
    sqrt2sigma_ = sqrt(2)*sigma;
    k_  = indus_k (sigma_, rc_);
    k1_ = indus_k1(k_, sigma);
    k2_ = indus_k2(k_, sigma_, rc_);
    return;
  }
  virtual double calculate(const Vec3<double>& x) const{
    double dist = getDistance(x, position_, box_size_);
    double deltaR = rmax_ - dist;
    if(dist > rmax_ + rc_) return 0.0;
    if(dist < rmax_ - rc_) return 1.0;
    //need to get distance accounting for pbc
    return (k1_*erf((deltaR)/sqrt2sigma_) - k2_*(deltaR) - 0.5)*heaviside(rc_ - fabs(deltaR)) + heaviside(rc_ + deltaR);
  }
  void updatePosition(const Vec3<double> pos){
    position_ = pos;
  }
private:
  Vec3<double> position_;
  double sigma_, sqrt2sigma_, rc_, rmax_, k_, k1_, k2_;
};

class cuboidalINDUSVolume : public INDUSVolume{
public:
  cuboidalINDUSVolume(const Vec3<double>& xmin, const Vec3<double>& xmax, double sigma, double cutoff){
    sigma_ = sigma;
    xmax_original_ = xmax;
    xmin_original_ = xmin;
    xc_ = cutoff;
    sqrt2sigma_ = sqrt(2)*sigma;
    k_  = indus_k (sigma_, xc_);
    k1_ = indus_k1(k_, sigma);
    k2_ = indus_k2(k_, sigma_, xc_);
    for(int i = 0; i < 3; i++){
      if(xmin_original_[i] < 0) xmin_original_[i] = 0;
    }
  }
  virtual void update_box_size(const Vec3<double>& box_size){
    box_size_ = box_size;
    //snap observation volume to edge of box
    for(int i = 0; i < 3; i++){
      if(xmin_original_[i] < 0) xmin_[i] = 0;
      else xmin_[i] = xmin_original_[i];
      if(xmax_original_[i] > box_size_[i]) xmax_[i] = box_size_[i];
      else xmax_[i] = xmax_original_[i];
    }
    return;
  }
  virtual double calculate(const Vec3<double>& x) const{
    //with the left boundary, needs closest pbc image to the left face
    //with the right boundary, needs closest pbc image to the right face
    //necessary step for observation volume near box edge
    Vec3<double> h_i = {0,0,0};
    for(int i = 0; i < 3; i++){
      double x_r, x_l, x_w; //right x, left x, wrapped x
      x_w = wrapNumber(x[i], box_size_[i]);
      if(x_w < xmin_[i] - xc_ || x_w > xmax_[i] + xc_) return 0.0;
      if(x_w > xmin_[i] + xc_ && x_w < xmax_[i] - xc_){
        h_i[i] = 1.0;
        continue;
      }
      x_r = getNearestImage1D(x[i], xmax_[i], box_size_[i]);
      x_l = getNearestImage1D(x[i], xmin_[i], box_size_[i]);
      double rbt = (k1_ * erf((xmax_[i]-x_r)/sqrt2sigma_) - k2_*(xmax_[i]-x_r) - 0.5)*heaviside(xc_ - fabs(xmax_[i]-x_r));
      double lbt = (k1_ * erf((x_l-xmin_[i])/sqrt2sigma_) - k2_*(x_l-xmin_[i]) - 0.5)*heaviside(xc_ - fabs(x_l-xmin_[i]));
      double hvt = heaviside(  xc_ + 0.5*(xmax_[i]-xmin_[i]) - fabs(x_w - 0.5*(xmin_[i] + xmax_[i]))  );
      h_i[i] = rbt + lbt + hvt;
    }
    return h_i[0] * h_i[1] * h_i[2];
  }
  protected:
  Vec3<double> position_;
  double sigma_, sqrt2sigma_, xc_, k_, k1_, k2_;
  Vec3<double> xmax_original_, xmin_original_, xmin_, xmax_;
};