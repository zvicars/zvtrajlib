#pragma once
#include <vector>
#include <limits>
#include <cmath>
static inline void makeHistogram(const std::vector<double>& data, double min_bin_, double max_bin_, double bin_size_, bool forceMin, bool forceMax, bool forceBS,
                                   std::vector<double>& x_vals, std::vector<int>& y_vals){
  if(!forceMin || !forceMax || !forceBS){
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    for(std::size_t i = 0; i < data.size(); i++){
        double datum = data[i];
        if(datum < min) min = datum;
        if(datum > max) max = datum;
    }
    min_bin_ = min;
    max_bin_ = max;
  }
  int nbins = std::round((max_bin_ - min_bin_) / bin_size_) + 1;
  //prepare ouput vectors
  x_vals.clear();
  x_vals.resize(nbins, 0.0);
  y_vals.clear();
  y_vals.resize(nbins, 0);
  for(int i = 0; i < nbins; i++){
      x_vals[i] = (double)i*bin_size_ + min_bin_;
  }

  for(std::size_t i = 0; i < data.size(); i++){
      double datum = data[i];
      int index = std::round((datum - min_bin_)/bin_size_);
      y_vals[index]++;
  }

  return;
}
static inline void makeHistogram2d(const std::vector<std::vector<double> >& data, double min_bin_, double max_bin_, double bin_size_, bool forceMin, bool forceMax, bool forceBS,
                                   std::vector<double> x_vals, std::vector<int> y_vals){
  //if I need to dynamically determine any bin bounds, loop through data to get min/max
  if(!forceMin || !forceMax || !forceBS){
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    for(std::size_t i = 0; i < data.size(); i++){
      for(std::size_t j = 0; j < data[i].size(); j++){
        double datum = data[i][j];
        if(datum < min) min = datum;
        if(datum > max) max = datum;
      }
    }
    if(!forceMin) min_bin_ = min;
    if(!forceMax) max_bin_ = max;
    if(!forceBS) bin_size_ = (max-min)/49.0; //default to a histogram with 50 points
  }

  int nbins = std::round((max_bin_ - min_bin_) / bin_size_) + 1;
  //prepare ouput vectors
  x_vals.clear();
  x_vals.resize(nbins, 0.0);
  y_vals.clear();
  y_vals.resize(nbins, 0);
  for(int i = 0; i < nbins; i++){
      x_vals[i] = (double)i*bin_size_ + min_bin_;
  }

  for(std::size_t i = 0; i < data.size(); i++){
    for(std::size_t j = 0; j < data[i].size(); j++){
      double datum = data[i][j];
      int index = std::round((datum - min_bin_)/bin_size_);
      y_vals[index]++;
    }
  }
  return;
}
