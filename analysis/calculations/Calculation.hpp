//A calculation takes in a full-frame's worth of information and calculates one or more things, for instance, N_tildeV and the corresponding forces. As calculations are persistent, they can also
//store the data from previous frames and output them at the end, though I may need to put more thought into it to prevent memory issues
//this is also probably where parallelizaton should be focused
#pragma once
#include "../../tools/Assert.hpp"
#include "../../tools/StringTools.hpp"
#include "../../interface/datatypes.hpp"
#include "../InputPack.hpp"
#include "../probevolumes/ProbeVolume.hpp"
#include "../atomgroups/AtomGroup.hpp"
#include <fstream>

class Calculation{
public:
  using KeyType = ParameterPack::KeyType;
  Calculation(InputPack& input);
  virtual ~Calculation(){
    return;
  }
  virtual void calculate(){return;}
  virtual std::string printConsoleReport(){return "";}
  virtual void finalOutput(){return;}
  virtual void update();
  virtual void output(){
    if(!doOutput()) return;
    printOutput();
    return;
  }
  bool isFinished(){
    if(box->time > end_){
      return 1;
    }
    return 0;
  }
  virtual void printOutput(){return;}
  virtual bool doOutput();
  virtual bool doCalculate();
  virtual bool doCalculateChild();
  bool hasUpdated(){
    return update_flag_;
  }
  bool hasCalculated(){
    return calculate_flag_;
  }
  void finish(){
    update_flag_ = 0;
    calculate_flag_ = 0;
  }
protected:
  std::string name_, type_, base_;
  double equilibration_, current_time_, end_;
  int output_freq_, calc_freq_, current_frame_; //number of frames to wait before outputting
  const Box* box = 0; //rather than passing in the box object each time, the input pack provides a pointer to this Calculation's box object
  bool update_flag_, calculate_flag_;
};