#include "AtomGroup_dynamic.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
AtomGroup_dynamic::AtomGroup_dynamic(InputPack& input) : AtomGroup{input} {
  input.params().readString("file", KeyType::Required, filename_);
  input.params().readFlag("permissive", KeyType::Required, allowOB_); 
  std::ifstream ifile(filename_);
  FANCY_ASSERT(ifile.is_open(), "Failed to open specified file " + filename_ + " for dynamic atomgroup " + name_ + ".");
  std::string line;
  while(std::getline(ifile, line)){
    if(line.length() == 0) continue;
    if(line.at(0) == '#') continue;
    std::stringstream ss(line);
    double time;
    int idx1;
    std::vector<int> indices;
    ss >> time;
    while(ss >> idx1){
      indices.push_back(idx1);
    }
    frames_.push_back(indices);
    times_.push_back(time);
  }
  FANCY_ASSERT(times_.size() != 0, "No indices provided to dynamic atomgroup " + name_ + ".");
  if(frames_.size() > 1){
    dt_ = times_[1] - times_[0];
    tmin_ = times_.front();
    tmax_ = times_.back();
  }
  return;
}
void AtomGroup_dynamic::update(){
  double time = box_->time;
  //if only one index is provided to avoid calling anything relying on dt, tmin, or tmax
  if(frames_.size() == 1){
    global_indices_ = frames_[0];
    return;
  }
  int index = round((time - tmin_) / dt_);
  if(index >= frames_.size()){
    if(!allowOB_){
      std::cout << "Asking for frames out of bounds of provided indices, \
      ending program. If you want to avoid this check, use the flag \"permissive = 1\" in the AtomGroup input file." << std::endl;
      exit(1);
    }
    global_indices_ = frames_.back();
  }
  else if(index < 0){
    if(!allowOB_){
      std::cout << "Asking for frames out of bounds of provided indices, \
      ending program. If you want to avoid this check, use the flag \"permissive = 1\" in the AtomGroup input file." << std::endl;
      exit(1);
    }
    global_indices_ = frames_.front();
  }
  else{
    global_indices_ = frames_[index];
  }
  return;
}