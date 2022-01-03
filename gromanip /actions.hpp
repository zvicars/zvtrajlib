#pragma once
#include "../interface/interface.hpp"
#include "GroManipData.hpp"
#include "boxtools.hpp"
#include <unordered_map>
//function primitive: void function(GroManipInputPack& input_data, const std::vector<std::string> command)
namespace boxtools::actions{
  std::unordered_map<std::string, void (*)(GroManipData&, const std::vector<std::string>&)> action_map;
  //use to initialize registry
  void registerActions(){
    return;
  }
  void merge(GroManipData&, const std::vector<std::string>&);
  void eraseres(GroManipData&, const std::vector<std::string>&);
  void duplicate(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyresname(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyatomname(GroManipData&, const std::vector<std::string>&);
  void rotate(GroManipData&, const std::vector<std::string>&);
  void flip(GroManipData&, const std::vector<std::string>&);
  void translate(GroManipData&, const std::vector<std::string>&);
  void loadgro(GroManipData&, const std::vector<std::string>&);
  void writegro(GroManipData&, const std::vector<std::string>&);
  void createprimitivevolume(GroManipData&, const std::vector<std::string>&);
  void createunionvolume(GroManipData&, const std::vector<std::string>&);
}