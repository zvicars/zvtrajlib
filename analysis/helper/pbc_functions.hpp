#pragma once
#include <array>
#include <vector>

static inline std::array<double,3> getNearestPeriodicPosition(const std::array<double,3>& pos, const std::array<double,3>& ref_pos, const std::array<double,3> box_size){
  std::array<double, 3> newpos, dx;
  newpos = pos;
  for(int i = 0; i < 3; i++){
    dx[i] = pos[i] - ref_pos[i];
    if(dx[i] > box_size[i]/2.0){
      newpos[i] -= box_size[i];
    } 
    else if(dx[i] < -box_size[i]/2.0){
      newpos[i] += box_size[i];
    } 
  }
  return newpos;
}