#pragma once
#include <map>
#include "../interface/datatypes.hpp"
#include "../tools/Assert.hpp"
#include <iostream>
class Volume;

class GroManipData{
public:
  GroManipData(){
    return;
  }
  ~GroManipData();
  Volume* findVolume(std::string name){
    auto it = vol_registry_.find(name);
    if(it != vol_registry_.end()) return it->second;
    return 0; //return nullptr if search fails
  }
  Box* findBox(std::string name){
    auto it = box_registry_.find(name);
    if(it != box_registry_.end()) return it->second;
    return 0; //return nullptr if search fails
  }
  void addVolume(std::string name, Volume* pv){
    vol_registry_.insert(std::pair<std::string, Volume*>{name, pv});
    return;
  }
  void addBox(std::string name, Box* box){
    box_registry_.insert(std::pair<std::string, Box*>{name, box});
    return;
  }
  void setVolumeRegistry(std::map<std::string, Volume*> vol_reg){
    vol_registry_ = vol_reg;
    return;
  }
  void setBoxRegistry(std::map<std::string, Box*> box_reg){
    box_registry_ = box_reg;
    return;
  }
  const std::map<std::string, Volume*>& VolumeMap() const {
    return vol_registry_; 
  }
  const std::map<std::string, Box*>& BoxMap() const {
    return box_registry_; 
  }

private:
  std::map<std::string, Volume*> vol_registry_; //contains pointer to the true map
  std::map<std::string, Box*> box_registry_; //contains pointer to the true map
};