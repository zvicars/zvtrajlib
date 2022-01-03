#pragma once
#include <map>
#include "../interface/datatypes.hpp"
#include "../tools/Assert.hpp"
class Volume;

class GroManipData{
public:
  GroManipData(){
    return;
  }
  ~GroManipData(){
    for(auto it = vol_registry_->begin(); it!=vol_registry_->end(); it++) {
      it = vol_registry_->erase(it);
    }
    for(auto it = box_registry_->begin(); it!=box_registry_->end(); it++) {
      it = box_registry_->erase(it);
    }
    delete vol_registry_;
    delete box_registry_;
  }
  Volume* findVolume(std::string name){
    auto it = vol_registry_->find(name);
    if(it != vol_registry_->end()) return it->second;
    return 0; //return nullptr if search fails
  }
  Box* findBox(std::string name){
    auto it = box_registry_->find(name);
    if(it != box_registry_->end()) return it->second;
    return 0; //return nullptr if search fails
  }
  void addVolume(std::string name, Volume* pv){
    if(vol_registry_ != 0) vol_registry_->insert(std::pair<std::string, Volume*>{name, pv});
    return;
  }
  void addBox(std::string name, Box* box){
    if(box_registry_ != 0) box_registry_->insert(std::pair<std::string, Box*>{name, box});
    return;
  }
  void initializeRegistries(){ //use this for master pack
    vol_registry_ = new std::map<std::string, Volume*>;
    box_registry_ = new std::map<std::string, Box*>;
  }
  void setVolumeRegistry(std::map<std::string, Volume*>* vol_reg){
    vol_registry_ = vol_reg;
    return;
  }
  void setBoxRegistry(std::map<std::string, Box*>* box_reg){
    box_registry_ = box_reg;
    return;
  }
  const std::map<std::string, Volume*>& VolumeMap() const {
    return *vol_registry_; 
  }
  const std::map<std::string, Box*>& BoxMap() const {
    return *box_registry_; 
  }

private:
  std::map<std::string, Volume*>* vol_registry_ = 0; //contains pointer to the true map
  std::map<std::string, Box*>* box_registry_ = 0; //contains pointer to the true map
};