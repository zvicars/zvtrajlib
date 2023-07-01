#pragma once
#include <cmath>
#include <array>
#include <complex>
#include <vector>
//implemented from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING
static inline double p_legendre(int l, int m, double x){
  int i,ll;
  double fact,oldfact,pll,pmm,pmmp1,omx2;
  if (m < 0 || m > l || abs(x) > 1.0) throw("Bad arguments in routine plgndr");
  pmm=1.0;
  if (m > 0) {
    omx2=(1.0-x)*(1.0+x);
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= omx2*fact/(fact+1.0);
      fact += 2.0;
    }
  }
  pmm=sqrt((2*m+1)*pmm/(4.0*M_PI));
  if (m & 1) pmm=-pmm;
  if (l == m) return pmm;
  else {
    pmmp1=x*sqrt(2.0*m+3.0)*pmm;
    if (l == (m+1)) return pmmp1;
    else { 
      oldfact=sqrt(2.0*m+3.0);
      for (ll=m+2;ll<=l;ll++) {
        fact=sqrt((4.0*ll*ll-1.0)/(ll*ll-m*m));
        pll=(x*pmmp1-pmm/oldfact)*fact;
        oldfact=fact;
        pmm=pmmp1;
        pmmp1=pll;
      }
      return pll;
    }
  }
  throw 0;
  return 0.0;
}
//non-normalized version
static inline double plgndr(int l, int m, double x){
  if (m < 0 || m > l || abs(x) > 1.0) throw("Bad args in routine plgndr");
  double prod = 1.0;
  for (int j=l-m+1;j<=l+m;j++) prod *= j;
  return sqrt(4.0*M_PI*prod/(2*l+1))*p_legendre(l,m,x);
}
static inline std::complex<double> Ylm(int l, int m, double costheta, double phi){
  double x = costheta;
  int m_eval = abs(m);
  std::complex<double> ylm_ret;
  ylm_ret = p_legendre(l, m_eval, x)*exp(std::complex<double>(0,1) * (double)m_eval * phi);
  if(m < 0) ylm_ret = std::pow(-1, m_eval)*std::conj(ylm_ret);
  return ylm_ret;
}
//use this for computing negatives of +m values
static inline std::complex<double> Ylm_invert(std::complex<double> ylm_in, double m){
  int m_pow = abs(m);
  return std::pow(-1, m_pow)*std::conj(ylm_in);
}