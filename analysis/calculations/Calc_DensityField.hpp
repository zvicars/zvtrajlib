#pragma once
#include "Calculation.hpp"
#include "../../tools/pbcfunctions.hpp"
class Calc_DensityField : public Calculation{
public:
  Calc_DensityField(InputPack& input);
  virtual void calculate();
  virtual void update();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
protected:
  Vec3<int> getIndex(const Vec3<double>& pos);
  bool getIndexNoWrap1D(double& pos, int dim, int& ret);
  int _map31(std::array<int,3>& coord) const{
    //pbc correct index
    for(int i = 0; i < 3; i++){
      coord[i] = wrapIndex(coord[i], npoints_[i]);
    }
    return (coord[2]*npoints_[0]*npoints_[1]) + (coord[1]*npoints_[0]) + coord[0];
  }
  Vec3<int> _map13(int idx) const{
    Vec3<int> eval;
    eval[2] = idx / (npoints_[0] * npoints_[1]);
    idx -= (eval[2] * npoints_[0] * npoints_[1]);
    eval[1] = idx / npoints_[0];
    eval[0] = idx % npoints_[0];
    return eval;
  }  
  bool isInBox(const Vec3<double>& x) const{
    if(x[0] < minx_[0] || x[0] > maxx_[0] || x[1] < minx_[1] || x[1] > maxx_[1] || x[2] < minx_[2] || x[2] > maxx_[2]) return 0;
    return 1;
  }
  double getGaussian(const Vec3<double>& position, const Vec3<int>& index);
  AtomGroup* atom_group_;
  std::vector<double> gridvals_, avggridvals_;
  Vec3<int> npoints_;
  Vec3<double> gridspacing_, avggridspacing_, minx_, maxx_, box_size_;
  bool hasBoxVec_, coarseGrain_;
  //stuff for coarse-graining
  double sigma_; Vec3<int> span_;
  int nframes_;
};