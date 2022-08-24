#include "tiles.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include "../tools/Matrix.hpp"
#include "../tools/smearfunctions.hpp"
#include "../tools/stlmath.hpp"
#include <sstream>
#include <iomanip>
double h_xy(Vec2<double> pos, Vec2<double> min, Vec2<double> max, double sigma, double cutoff){
  double eval = 1.0;
  for(int i = 0; i < 2; i++){
    eval *= h_x(pos[i], min[i], max[i], sigma, cutoff);
  }
  return eval;
}

int main(){
  Tile2d ht;
  double s = 6.4;
  ht.positions.push_back({0,0});
  ht.positions.push_back({0.5, 1.0/6.0});
  ht.positions.push_back({0.5, 0.5});
  ht.positions.push_back({0, 2.0/3.0});
  ht.spacing = {sqrt(3.0)*s, 3.0*s};

  //need to fill a 50x50 grid
  Vec2<int> ntiles = {ceil(50.0 / (sqrt(3.0)*s))  + 1 , ceil(50.0 / (3.0*s)) + 1};
  std::vector<Vec2<double> > positions_final;
  for(int i = -1; i < ntiles[0]; i++){
    for(int j = -1; j < ntiles[1]; j++){
      auto positions = ht.generate_offset_positions({i,j});
      positions_final.insert(positions_final.end(), positions.begin(), positions.end());
    }
  }

  //generate gridpoints
  //software takes in a 50x50 image normalized from 0 to 1
  for(double sigma_ = 1.5; sigma_ < 8.0; sigma_+=0.1){
    Matrix2d<double, 50, 50> grid_;
    grid_.fill(0.0);
    double gs_ = 1.0;
    //double sigma_ = 2.0;
    double cutoff_ = 3.0*sigma_;
    for(auto pos : positions_final){
      int imin = std::max((int)floor(pos[0] - cutoff_ - 1), 0);
      int imax = std::min((int)floor(pos[0] + cutoff_ + 1), 50);
      int jmin = std::max((int)floor(pos[1] - cutoff_ - 1), 0);
      int jmax = std::min((int)floor(pos[1] + cutoff_ + 1), 50);    
      for(int i = imin; i < imax; i++){
        for(int j = jmin; j < jmax; j++){
          Vec2<double> min = {i*gs_, j*gs_};
          Vec2<double> max = {i*gs_ + gs_, j*gs_ + gs_};
          Vec2<int> idx = {i,j};
          grid_.at(idx) += h_xy(pos, min, max, sigma_, cutoff_);
        }
      }    
    }
    double max = -1.0;
    double scale = 1.0;
    for(int i = 0; i < grid_.size_1d(); i++){
      if(grid_.read_at_1d(i) > max) max = grid_.read_at_1d(i);
    }
    //set the peak value to a given number
    for(int i = 0; i < grid_.size_1d(); i++){
      grid_.at_1d(i) *= scale*0.84/0.0404505;
    }  
    std::stringstream oss;
    oss << std::setprecision(3) << sigma_ << ".txt";
    std::cout << "\'" << std::setprecision(3) << sigma_ << "\', ";
    std::ofstream ofile(oss.str());
    for(int i = 0; i < 50; i++){
      for(int j = 0; j < 50; j++){
        Vec2<int> idx = {i,j};
        ofile << grid_.read_at(idx) << " ";
      }
      ofile << "\n";
    }
    ofile.close();
  }
  return 0;
}