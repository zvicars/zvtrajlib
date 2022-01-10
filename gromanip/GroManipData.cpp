 #include "GroManipData.hpp"
 #include "volume.hpp"
GroManipData::~GroManipData()
{
  for(auto it = vol_registry_.begin(); it!=vol_registry_.end(); it++) {
    auto ptr = it->second;
    delete ptr;
  }
  for(auto it = box_registry_.begin(); it!=box_registry_.end(); it++) {
    auto ptr = it->second;
    delete ptr;
  }
  return;
}