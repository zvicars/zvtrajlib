//Zachariah Vicars
//8/6/2022
//Computing H-bond profile for systems
//Reads in trajectory and outputs all of the h-bonded waters
//h-bonds are defined as O-H--O triplets with -- < r_h and an an angle of less than r_cut
//stores results as a vector of h-bond objects for each frame that can be referenced for other calculations
#pragma once
#include "Calculation.hpp"
#include <unordered_set>
struct HBond{
  Vec3<double> x, v; //origin and direction of the hydrogen bond
  double r, angle; //computed angle
  std::array<int, 2> indices; //O and H index
  std::array<double,2> charges;
};
class Calc_HBonds : public Calculation{
  public:
    Calc_HBonds(InputPack& input);
    virtual void calculate();
    virtual void finalOutput();
    virtual void update();
    virtual void output();
    const std::vector<HBond>& getHBonds() const{
      return hbonds_;
    }
  protected:
    AtomGroup* atom_group_, * oxygen_group_;
    double r_, angle_, r_oh_;
    std::vector<HBond> hbonds_;
    Vec3<double> box_size_;
    //valid formats: tip4p, spc/e
    int natpermol_;
    Vec<int> rel_h_indices_;
    int rel_charge_index_;
    bool write_traj_;
};