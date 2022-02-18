#include "PV_DynamicBox.hpp"
PV_DynBox::PV_DynBox(InputPack& input):ProbeVolume{input}
{
  Vec<double> x_range, y_range, z_range;
  axis_ranges_.resize(6, 0.0);
  input.params().readVector("x_range", KeyType::Required, x_range);
  input.params().readVector("y_range", KeyType::Required, y_range);
  input.params().readVector("z_range", KeyType::Required, z_range);
  box = input.getBox();
  FANCY_ASSERT(x_range.size() == 2, "Invalid x-range given for DynBox PV.");
  FANCY_ASSERT(y_range.size() == 2, "Invalid y-range given for DynBox PV.");
  FANCY_ASSERT(z_range.size() == 2, "Invalid z-range given for DynBox PV.");
  for(int i = 0; i < 2; i++){
    axis_ranges_[i*3] = x_range[i];
    axis_ranges_[i*3+1] = y_range[i];
    axis_ranges_[i*3+2] = z_range[i];
  }
  return;
}

double PV_DynBox::compute(Vec3<double> position){
  for(int i = 0; i < 3; i++){
    if( position[i] < axis_ranges_[i] + xcom_[i] || position[i] > axis_ranges_[i+3] + xcom_[i] ) return 0.0; 
  }
  return 1.0;
}

void PV_DynBox::update(){
  //place the probe volume at the center of mass of the provided atomgroup
  int natoms = atomgroup_->getIndexCount();
  xcom_ = {0.0};
  for(int i = 0; i < natoms; i++){
    for(int i = 0; i < 3; i++){
      xcom_[i] += box_->getAtomPosition(atomgroup_->getIndex(i))[i];
    }
  }
  for(int i = 0; i < 3; i++){
    xcom_[i] *= 1.0/(double)natoms;
  }
  update_flag_ = 1;
  return;
}