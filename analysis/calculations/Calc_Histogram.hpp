#pragma once
#include "Calculation.hpp"
#include "../helper/make_histogram.hpp"
class Calculation_Histogram : public Calculation{
public:
  Calculation_Histogram(InputPack& input);
  virtual ~Calculation_Histogram(){
    return;
  }
  virtual void calculate(){
    Calculation::calculate();
    return;
  }
  virtual std::string printConsoleReport(){
    return "";
  }
  virtual void finalOutput(){
    Calculation::finalOutput();
    return;
  }
  virtual void update(){
    Calculation::update();
    return;
  };
  virtual void output(){
    Calculation::output();
    return;
  }
protected:
  double bin_size_, min_bin_, max_bin_; //standard histogram output options, min/max bin can be found dynamically too
  bool doHistogram=0, doTimeseries=0, forceMin=0, forceMax=0, forceBS=0;
};