#pragma once
#include <vector>
#include <array>
#include <cmath>
static inline double wrapNumber(double x, double x_max){
  return x - x_max * floor( x / x_max );
}
static inline int wrapIndex(int i, int i_max) {
   return ((i % i_max) + i_max) % i_max;
}
//pbc implementation using floor function that wraps all numbers to box 
static inline void placeInsideBox(std::array<double, 3>& position, const std::array<double,3>& box_size){
  for(int i = 0; i < 3; i++){
    position[i] = position[i] - box_size[i]*floor(position[i] / box_size[i]);
  }
  return;
}
static inline void placeInsideBox(std::array<double, 3>& position, const std::array<double,9>& box_size){
  //same function as before, but works for non-orthorhombic boxes
  throw 0;
  return;
}
static inline double getDistanceNoPBC(const std::array<double, 3>& x1, const std::array<double, 3>& x2){
  double dist = 0.0;
  for(int i = 0; i < 3; i++){
    dist += (x1[i]-x2[i])*(x1[i]-x2[i]);
  }
  return std::sqrt(dist);
}
//handles arbitrary number of pbc wraps
static inline double getDistance(const std::array<double, 3>& x1, const std::array<double, 3>& x2, const std::array<double,3>& box_size){
  std::array<double,3> dx;
  for(int i = 0; i < 3; i++){
    dx[i] = x2[i] - x1[i];
    while(dx[i] >= 0.5*box_size[i]){
      dx[i] -= box_size[i];
    }
    while(dx[i] < -0.5*box_size[i]){
      dx[i] += box_size[i];
    }    
  }
  return std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
}
static inline double getDistance(const std::array<double, 3>& x1, const std::array<double, 3>& x2, const std::array<double,9>& box_size){
  throw 0;
  return 0.0;
}
static inline double getNearestImage1D(double x, double xref, double box_size){
  while(x - xref >= 0.5*box_size){
    x -= box_size;
  }
  while(x - xref < -0.5*box_size){
    x += box_size;
  }
  return x; 
}

static inline void getNearestImage3D(std::array<double,3>& x, const std::array<double,3>& xref, const std::array<double,3>& box_size){
  for(int i = 0; i < 3; i++){
    x[i] = getNearestImage1D(x[i], xref[i], box_size[i]);
  }
  return; 
}