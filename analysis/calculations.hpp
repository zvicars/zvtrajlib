//A calculation takes in a full-frame's worth of information and calculates one or more things, for instance, N_tildeV and the corresponding forces. As calculations are persistent, they can also
//store the data from previous frames and output them at the end, though I may need to put more thought into it to prevent memory issues
//this is also probably where parallelizaton should be focused
#pragma once
#include "../io/parser.hpp" //needs to be able to register itself
#include "../interface/datatypes.hpp"

class Calculation{
public:
  Calculation(std::string input){
    viaKey<std::string>("name", input, name_);
    return;
  }
  virtual void calculate(const Box& box) = 0;
  virtual void write(std::string base) = 0;
private:
  std::string name_;
};


class Calc_Nv{
public:
  Calc_Nv(std::string input){
    std::string pv_name_;
    viaKey<std::string>("probe_volume", input, pv_name_);
    pv_ = reg::probe_volumes.findObject("pv_name_");
  }
  void calculate(const Box& box){
    for(int i = 0; i < box.atoms.size(); i++){
      pv_->calculate(box.atoms[i].x);
    }
  }
  
private:
  Vec<double> time;
  Vec<int> step;
  Vec<int> count;
  std::string pv_name_;
  ProbeVolume* pv_;
};