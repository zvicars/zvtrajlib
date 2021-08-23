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
  virtual void calculate()=0;
  virtual std::string printConsoleReport()=0;
  virtual void finalOutput(){return;}
  virtual void update();
  virtual void output(){
    if(!doOutput()) return;
    printOutput();
    return;
  }
protected:
  virtual void printOutput(){return;}
  //check to see if this is a step where you should output
  virtual bool doOutput();
  virtual bool doCalculate();
  std::string name_, type_, base_;
  double equilibration_, current_time_;
  int output_freq_, calc_freq_, current_frame_; //number of frames to wait before outputting
  const Box* box = 0; //rather than passing in the box object each time, the input pack provides a pointer to this Calculation's box object
  //options for final output modes, typical behavior is just to spit out means/variances, but time series data and histogram data can be toggled
  //here if the calculation allows it all of this is handled in the "final output" as it is expected that this data will be sufficiently small
  //to be stored in memory and outputted at the end
  double bin_size_, min_bin_, max_bin_; //standard histogram output options, min/max bin can be found dynamically too
  bool doHistogram=0, doTimeseries=0, forceMin=0, forceMax=0, forceBS=0;
};