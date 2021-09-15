//calculate end group dipole for looking at dipole-surface interactions
#pragma once
#include "Calculation.hpp"
#include <cmath>
class Calc_Dipole : public Calculation{
public:
  Calc_Dipole(InputPack& input);
  virtual void calculate();
  virtual void update();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:
  AtomGroup* atom_group_;
  //expecting a sequential list of molecules, will modulo with molsize to get appropriate charge for atom
  int natoms_, nmols_, at_per_mol_;
  Vec<double> masses_, charges_;
  Vec3<double> normal_; //normal vector of surface
  Vec<Vec3<double> > dipoles_;
  Vec<double> costhetas_, times_, steps_;
  Vec<bool> included_atoms_;
  double costheta_avg_frame_, costheta_avg_total_, var_;
  int costheta_avg_total_counter_;
};