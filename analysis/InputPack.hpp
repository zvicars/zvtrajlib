#pragma once
#include <map>
#include "../tools/InputParser.hpp"
#include "../tools/Stacktrace.hpp"
#include "../interface/datatypes.hpp"
class ProbeVolume;
class Calculation;
class AtomGroup;

//input packs contain a parameterpack and pointer to the registries that different calculations might need to functon,
class InputPack{
public:
  InputPack(){
    return;
  }
  InputPack(const ParameterPack* params, const Box* box){
    params_ = params;
    initializeRegistries();
    box_ = box;
    return;
  }
  ~InputPack(){
    if(is_master_pack){ //only a pack that initialized the registries should delete them
      for(auto it = pv_registry_->begin(); it!=pv_registry_->end(); it++) {
        it = pv_registry_->erase(it);
      }
      for(auto it = calc_registry_->begin(); it!=calc_registry_->end(); it++) {
       it = calc_registry_->erase(it);
      }
      for(auto it = ag_registry_->begin(); it!=ag_registry_->end(); it++) {
       it = ag_registry_->erase(it);
      }
      delete pv_registry_;
      delete calc_registry_;
      delete ag_registry_;
    }
  }
  const Box* getBox(){
    return box_;
  }
  ProbeVolume* findProbeVolume(std::string name){
    auto it = pv_registry_->find(name);
    if(it != pv_registry_->end()) return it->second;
    return 0; //return nullptr if search fails
  }

  Calculation* findCalculation(std::string name){
    auto it = calc_registry_->find(name);
    if(it != calc_registry_->end()) return it->second;
    return 0; //return nullptr if search fails
  }
  AtomGroup* findAtomGroup(std::string name){
    auto it = ag_registry_->find(name);
    if(it != ag_registry_->end()) return it->second;
    return 0; //return nullptr if search fails
  }
  void addProbeVolume(std::string name, ProbeVolume* pv){
    if(pv_registry_ != 0) pv_registry_->insert(std::pair<std::string, ProbeVolume*>{name, pv});
    return;
  }

  void addCalculation(std::string name, Calculation* calc){
    if(calc_registry_ != 0) calc_registry_->insert(std::pair<std::string, Calculation*>{name, calc});
    return;
  }
  void addAtomGroup(std::string name, AtomGroup* ag){
    if(ag_registry_ != 0) ag_registry_->insert(std::pair<std::string, AtomGroup*>{name, ag});
    return;
  }
  void initializeRegistries(){ //use this for master pack
    pv_registry_ = new std::map<std::string, ProbeVolume*>;
    calc_registry_ = new std::map<std::string, Calculation*>;
    ag_registry_ = new std::map<std::string, AtomGroup*>;
  }

  void setProbeVolumeRegistry(std::map<std::string, ProbeVolume*>* pv_reg){
    pv_registry_ = pv_reg;
    return;
  }
  void setCalculationRegistry(std::map<std::string, Calculation*>* calc_reg){
    calc_registry_ = calc_reg;
    return;
  }
  void setAtomGroupRegistry(std::map<std::string, AtomGroup*>* ag_reg){
    ag_registry_ = ag_reg;
    return;
  }
  const ParameterPack& params() const {
    return *params_; 
  }
  void setParams(const ParameterPack* param){
    params_ = param;
  }
  void setBox(const Box* box){
    box_ = box;
  }
  const std::map<std::string, ProbeVolume*>& ProbeVolumeMap() const {
    return *pv_registry_; 
  }
  const std::map<std::string, Calculation*>& CalculationMap() const {
    return *calc_registry_; 
  }  

  const std::map<std::string, AtomGroup*>& AtomGroupMap() const {
    return *ag_registry_; 
  }  
  std::vector<InputPack> buildDerivedInputPacks(std::string key){
    std::vector<InputPack> inputpacks;
    auto parameterpacks = params().findParameterPacks(key, ParameterPack::KeyType::Optional);
    inputpacks.resize(parameterpacks.size());
    for(std::size_t i = 0; i < parameterpacks.size(); i++){
      inputpacks[i].setCalculationRegistry(calc_registry_);
      inputpacks[i].setProbeVolumeRegistry(pv_registry_);
      inputpacks[i].setAtomGroupRegistry(ag_registry_);
      inputpacks[i].setParams(parameterpacks[i]);
      inputpacks[i].setBox(box_);
    }
    return inputpacks;
  }


private:
  const ParameterPack* params_;
  std::map<std::string, ProbeVolume*>* pv_registry_ = 0; //contains pointer to the true map
  std::map<std::string, Calculation*>* calc_registry_ = 0; //contains pointer to the true map
  std::map<std::string, AtomGroup*>* ag_registry_ = 0; //contains pointer to the true map
  const Box* box_;
  bool is_master_pack = 0;
};