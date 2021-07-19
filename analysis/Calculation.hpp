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

class Calculation{
public:
  using KeyType = ParameterPack::KeyType;
  Calculation(InputPack& input);
  virtual ~Calculation(){
    return;
  }
  virtual void calculate(const Box& box)=0;
  virtual std::string printConsoleReport()=0;
  //virtual void write(std::string base) = 0;
protected:
std::string name_, type_;
private:
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

