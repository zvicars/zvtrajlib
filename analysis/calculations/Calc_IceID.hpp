//uses a variety of criteria to compute an approximate number of true ice-like molecules within a probe volume based on one or more criteria
//this is then output as an xtc file on a per-step basis and computes statistics on it

#include "Calc_Histogram.hpp"
#include "../../interface/interface.hpp"
class Calc_IceID : public Calculation_Histogram{
public:
  Calc_IceID(InputPack& input);
  virtual void update();
  virtual void calculate();
  virtual void output();
  virtual void finalOutput();
protected:
  void performIteration(std::vector<int>& ice_indices);
  //atomgroups and probevolumes
  AtomGroup* ice_group_, * surface_group_, * water_group_;
  std::string pv_name_, icename_, surfname_, watername_;
  ProbeVolume* pv_;
  //parameters
  double ice_radius_, surface_radius_;
  int ice_thresh_, surf_thresh_;
  std::vector<int> final_ice_indices_;
  int n_, nrep_;
  std::vector<double> n_vec_;
  std::vector<int> step_vec_;
  std::vector<double> t_vec_;
  //xtc output
  xdr::XDRFILE* output_handle_;
  int xdr_natoms_, xdr_step_;
  float xdr_time_;
  xdr::matrix xdr_box_;
  xdr::rvec* xdr_x_;
  float xdr_prec_;
  bool initialized_;
  std::string filename_;
};