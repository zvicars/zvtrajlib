//A calculation takes in a full-frame's worth of information and calculates one or more things, for instance, N_tildeV and the corresponding forces. As calculations are persistent, they can also
//store the data from previous frames and output them at the end, though I may need to put more thought into it to prevent memory issues
//this is also probably where parallelizaton should be focused
#pragma once
#include "../tools/Assert.hpp"
#include "../tools/GenericFactory.hpp"
#include "../interface/datatypes.hpp"
#include "InputPack.hpp"
#include "ProbeVolume.hpp"


class Calculation{
public:
  using KeyType = ParameterPack::KeyType;
  Calculation(InputPack& input);
  virtual void calculate(const Box& box) = 0;
  virtual void write(std::string base) = 0;
private:
  std::string name_;
};

class Calc_Nv : public Calculation{
public:
  Calc_Nv(InputPack& input);
  void calculate(const Box& box);
private:
  Vec<double> time;
  Vec<int> step;
  Vec<int> count;
  std::string pv_name_;
  ProbeVolume* pv_;
};


namespace CalculationRegistry {

using Key  = std::string;
using Base = Calculation;
using Factory = GenericFactory<Base, Key, InputPack&>;

template<typename D>
using Register = RegisterInFactory<Base, D, Key, InputPack&>;
}