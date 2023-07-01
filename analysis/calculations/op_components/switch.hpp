#pragma once
#include "../../../tools/Assert.hpp"
#include "../../../tools/stlmath.hpp"
#include "../../../interface/datatypes.hpp"
#include "../../../tools/smearfunctions.hpp"
#include "../../../tools/pbcfunctions.hpp"
//basic switching function class using INDUS-style h(x) one-boundary expression
class simpleSwitch{
public:
  simpleSwitch(double xmax, double sigma, double cutoff){
    xmax_ = xmax;
    sigma_ = sigma;
    xc_ = cutoff;
    sqrt2sigma_ = sqrt(2)*sigma;
    k_  = indus_k (sigma_, xc_);
    k1_ = indus_k1(k_, sigma);
    k2_ = indus_k2(k_, sigma_, xc_);
    return;
  }
  double calculate(double x) const{
    if(x < xmax_ - xc_) return 0.0;
    if(x > xmax_ + xc_) return 1.0;
    //need to get distance accounting for pbc
    return (k1_*erf((x-xmax_)/sqrt2sigma_) - k2_*(x-xmax_) - 0.5)*heaviside(xc_ - fabs(x-xmax_)) + heaviside(xc_ + x - xmax_);
  }
private:
  double sigma_, sqrt2sigma_, xc_, xmax_, k_, k1_, k2_;
};

class simpleTwoSidedSwitch{
public:
  simpleTwoSidedSwitch(double xmin, double xmax, double sigma, double cutoff){
    sigma_ = sigma;
    xmax_ = xmax;
    xmin_ = xmin;
    xc_ = cutoff;
    sqrt2sigma_ = sqrt(2)*sigma;
    k_  = indus_k (sigma_, xc_);
    k1_ = indus_k1(k_, sigma);
    k2_ = indus_k2(k_, sigma_, xc_);
  }
  double calculate(double x) const{
    if(x < xmin_ - xc_ || x > xmax_ + xc_) return 0.0;
    if(x > xmin_ + xc_ && x < xmax_ - xc_) return 1.0;
    double rbt = (k1_ * erf((xmax_-x)/sqrt2sigma_) - k2_*(xmax_-x) - 0.5)*heaviside(xc_ - fabs(xmax_-x));
    double lbt = (k1_ * erf((x-xmin_)/sqrt2sigma_) - k2_*(x-xmin_) - 0.5)*heaviside(xc_ - fabs(x-xmin_));
    double hvt = heaviside(  xc_ + 0.5*(xmax_-xmin_) - fabs(x - 0.5*(xmin_ + xmax_))  );
    return rbt + lbt + hvt;
  }
private:
  double xmax_, xmin_;
  double sigma_, sqrt2sigma_, xc_, k_, k1_, k2_;
};