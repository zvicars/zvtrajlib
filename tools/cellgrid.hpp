//simple algorithm for generating an nxnxn cell_grid
#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include "pbcfunctions.hpp"
//partitions the simulation cell into different regions where neigbors within cutoff can only be in neigboring grid cells
class CellGrid{
public:
  CellGrid(){
    return;
  }
  CellGrid(double cutoff, std::array<double, 3> box_size){
    for(int i = 0; i < 3; i++){
      nx[i] = ceil(box_size[i]/cutoff);
      sxi[i] = (double)nx[i]/box_size[i];
    }
    positions.resize(nx[0]*nx[1]*nx[2]);
    return;
  }
  void reset(double cutoff, std::array<double, 3> box_size){
    positions.clear();
    for(int i = 0; i < 3; i++){
      nx[i] = ceil(box_size[i]/cutoff);
      sxi[i] = (double)nx[i]/box_size[i];
    }
    positions.resize(nx[0]*nx[1]*nx[2]);
    return;    
  }
  std::vector<int> getNearbyIndices(const std::array<double, 3>& position) const{
    //searches 27 nearest cells and returns a list of indices stored within those cells
    std::vector<int> output;
    std::array<int,3> coord, new_coord;
    for(int i = 0; i < 3; i++){
      coord[i] = position[i]*sxi[i];
    }
    for(int i = -1; i <= 1; i++){
      new_coord[0] = coord[0] + i;
      for(int j = -1; j <= 1; j++){
        new_coord[1] = coord[1] + j;
        for(int k = -1; k <= 1; k++){
          new_coord[2] = coord[2] + k;
          int idx1d = _map31(new_coord);
          output.insert(output.end(), positions[idx1d].begin(), positions[idx1d].end());
        }
      }
    }
    return output;
  }
  //removes a reference index from the list to prevent an object from checking a distance with itself
  std::vector<int> getNearbyIndices(int refIdx, const std::array<double, 3>& position) const{
    //searches 27 nearest cells and returns a list of indices stored within those cells
    std::vector<int> output;
    std::array<int,3> coord, new_coord;
    for(int i = 0; i < 3; i++){
      coord[i] = position[i]*sxi[i];
    }
    for(int i = -1; i <= 1; i++){
      new_coord[0] = coord[0] + i;
      for(int j = -1; j <= 1; j++){
        new_coord[1] = coord[1] + j;
        for(int k = -1; k <= 1; k++){
          new_coord[2] = coord[2] + k;
          int idx1d = _map31(new_coord);
          if(positions[idx1d].size() > 0) output.insert(output.end(), positions[idx1d].begin(), positions[idx1d].end());
        }
      }
    }
    output.erase(remove(output.begin(), output.end(), refIdx), output.end());
    return output;
  }
  void addIndexToGrid(int atom_index, const std::array<double, 3>& position){
    std::array<int,3> coord;
    for(int i = 0; i < 3; i++){
      coord[i] = position[i]*sxi[i];
    }
    int index = _map31(coord);
    positions[index].push_back(atom_index);
    return;
  }
private:
  int _map31(std::array<int,3>& coord) const{
    //pbc correct index
    for(int i = 0; i < 3; i++){
      coord[i] = wrapIndex(coord[i], nx[i]);
    }
    return (coord[2]*nx[0]*nx[1]) + (coord[1]*nx[0]) + coord[0];
  }
  std::vector<std::vector<int> > positions; //nx*ny*nz
  std::array<int,3> nx;
  std::array<double,3> sxi; //inverse of cell sizes
};