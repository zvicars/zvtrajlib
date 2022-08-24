//simple algorithm for generating an nxnxn cell_grid
#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include "pbcfunctions.hpp"
#include "Matrix.hpp"
#include <pthread.h>
//partitions the simulation cell into different regions where neigbors within cutoff can only be in neigboring grid cells
class CellGrid{
public:
  CellGrid(){
    return;
  }
  CellGrid(double cutoff, std::array<double, 3> box_size, bool pbc = 1){
    for(int i = 0; i < 3; i++){
      nx[i] = ceil(box_size[i]/cutoff);
      sxi[i] = (double)nx[i]/box_size[i];
    }
    positions.initialize(nx);
    mutex_matrix.initialize(nx); 
    pbc_ = pbc;
    return;
  }
  void reset(double cutoff, std::array<double, 3> box_size){
    for(int i = 0; i < 3; i++){
      nx[i] = ceil(box_size[i]/cutoff);
      sxi[i] = (double)nx[i]/box_size[i];
    }
    positions.reset(nx);
    mutex_matrix.reset(nx);
    return;    
  }
  const std::vector<int>& getNearbyIndices(const std::array<double, 3>& position) const{
    //searches 27 nearest cells and returns a list of indices stored within those cells
    std::vector<int> output;
    std::array<int,3> coord, new_coord;
    for(int i = 0; i < 3; i++){
      coord[i] = wrapIndex(floor(position[i]*sxi[i]), nx[i]);
    }
    return positions.read_at(coord);
  }
  //removes a reference index from the list to prevent an object from checking a distance with itself
  std::vector<int> getNearbyIndices(int refIdx, const std::array<double, 3>& position) const{
    //searches 27 nearest cells and returns a list of indices stored within those cells
    std::vector<int> output;
    std::array<int,3> coord, new_coord;
    for(int i = 0; i < 3; i++){
      coord[i] = wrapIndex(floor(position[i]*sxi[i]), nx[i]);
    }
    output = positions.read_at(coord); 
    output.erase(remove(output.begin(), output.end(), refIdx), output.end());
    return output;
  }
  void addIndexToGrid(int atom_index, const std::array<double, 3>& position){
    std::array<int,3> coord;
    for(int i = 0; i < 3; i++){
      coord[i] = floor(position[i]*sxi[i]);
      if(coord[i] < 0 || coord[i] >= nx[i]){
        if(!pbc_) return;
        else{
          coord[i] = wrapIndex(coord[i], nx[i]);
        }
      }
    }
    for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++)
    for(int k = -1; k <= 1; k++){
      std::array<int,3> new_coord = {coord[0] + i, coord[1] + j, coord[2] + k};
      bool isOB = 0;
      for(int w = 0; w < 3; w++){
        if(new_coord[w] < 0 || new_coord[w] >= nx[w]){
          if(pbc_) new_coord[w] = wrapIndex(new_coord[w], nx[w]);
          else isOB = 1;
        }
      }
      if(isOB) continue;
      pthread_mutex_t* mt = &mutex_matrix.at(new_coord);
      pthread_mutex_lock(mt);
      positions.at(new_coord).push_back(atom_index);
      pthread_mutex_unlock(mt);
    }

    return;
  }
private:
  Matrix<std::vector<int>, 3> positions;
  Matrix<pthread_mutex_t, 3> mutex_matrix; 
  std::array<std::size_t,3> nx;
  std::array<double,3> sxi; //inverse of cell sizes
  bool pbc_ = 1;
};