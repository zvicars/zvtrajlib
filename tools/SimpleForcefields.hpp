//some simple forcefields, LJ, Couloumb interactions, etc.
//Zachariah Vicars, 7-7-2022
#pragma once
#include <cmath>
static inline double LJ_6_12(double r, double epsilon, double sigma){
  double sigma6 = std::pow(sigma/r, 6.0);
  return -4.0*epsilon*(sigma6 - sigma6*sigma6);
}