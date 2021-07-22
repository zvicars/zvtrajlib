//A calculation takes in a full-frame's worth of information and calculates one or more things, for instance, N_tildeV and the corresponding forces. As calculations are persistent, they can also
//store the data from previous frames and output them at the end, though I may need to put more thought into it to prevent memory issues
//this is also probably where parallelizaton should be focused
#pragma once
#include "../tools/Assert.hpp"
#include "../tools/GenericFactory.hpp"
#include "../tools/StringTools.hpp"
#include "../interface/datatypes.hpp"
#include "InputPack.hpp"
#include "ProbeVolume.hpp"
#include "AtomGroup.hpp"
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
  virtual void printOutput(){return;}
  virtual void finalOutput(){return;}
  virtual bool doCalculate();
  virtual bool doOutput();
protected:
  std::string name_, type_, base_;
  double equilibration_, current_time_;
  int output_freq_, calc_freq_, current_frame_; //number of frames to wait before outputting
  const Box* box = 0; //rather than passing in the box object each time, the input pack provides a pointer to this Calculation's box object
  //options for final output modes, typical behavior is just to spit out means/variances, but time series data and histogram data can be toggled
  //here if the calculation allows it all of this is handled in the "final output" as it is expected that this data will be sufficiently small
  //to be stored in memory and outputted at the end
  double bin_size_, min_bin_, max_bin_; //standard histogram output options, min/max bin can be found dynamically too
  bool doHistogram=0, doTimeseries=0, forceMin=0, forceMax=0, forceBS;
};

namespace CalculationRegistry {

// Manages the creation of Steinhardt objects
using Factory = GenericFactory<
  Calculation,   // base class
  std::string,  // key type
  InputPack&  // input types
>;

// Object that registers the mapping with the factory
template<typename S>
using Register = RegisterInFactory<
  Calculation,   // base class
  S,            // derived class
  std::string,  // key type
  InputPack&
>;
} 

