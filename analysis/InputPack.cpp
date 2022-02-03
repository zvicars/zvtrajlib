#include "InputPack.hpp"
#include "probevolumes/ProbeVolume.hpp"
#include "atomgroups/AtomGroup.hpp"
#include "calculations/Calculation.hpp"
InputPack::~InputPack()
{
  if(is_master_pack){ //only a pack that initialized the registries should delete them
    for(auto it = pv_registry_->begin(); it!=pv_registry_->end();) {
      delete it->second;
      it = pv_registry_->erase(it);
    }
    for(auto it = calc_registry_->begin(); it!=calc_registry_->end();) {
      delete it->second;
      it = calc_registry_->erase(it);
    }
    for(auto it = ag_registry_->begin(); it!=ag_registry_->end();) {
      delete it->second;
      it = ag_registry_->erase(it);
    }
    delete pv_registry_;
    delete calc_registry_;
    delete ag_registry_;
  }
}
