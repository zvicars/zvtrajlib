#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
//attempt to calculate histogram in position/momenta space for each voxel in the grid
class Calc_EntropyProfile : public Calc_DensityField{
  public:
    Calc_EntropyProfile(InputPack& input):Calc_DensityField{input}{
      point_data_.resize(gridvals_.size()); 
      return;
    }
    virtual void update(){
      Calc_DensityField::update();
      if(firstUpdate){
        FANCY_ASSERT(box->hasMasses && box->hasVelocities,
        "lacking sufficient information to calculate entropy profile. Need particle masses and topologies");
      }
      return;
    }
    virtual void calculate(){
      if(!doCalculate()) return;
      Calc_DensityField::calculate();
      auto indices = atom_group_->getIndices();
      for(auto index : indices){
          auto& atom = box->atoms[index];
          auto position = atom.x;
          placeInsideBox(position, box_size_);
          auto idx_ref = getIndex(position);
          std::array<double,6> temp = {position[0], position[1], position[2],
          atom.v[0]*atom.mass, atom.v[1]*atom.mass, atom.v[2]*atom.mass};
          point_data_[_map31(idx_ref)].push_back(temp);
      }
        return;
    }
    virtual void finalOutput(){
      Calc_DensityField::finalOutput();
      return;
    }
  protected:
  //x, y, z, px, py, pz
  //using maximimum likelihood arguments to get best fit for entropy of a gridcell
  std::vector< std::vector< std::array<double, 6> > > point_data_;
  bool firstUpdate = 1;
};