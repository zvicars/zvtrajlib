#pragma once
#include "../interface/interface.hpp"
#include "GroManipData.hpp"
#include "boxtools.hpp"
#include <unordered_map>
#include <iostream>
//function primitive: void function(GroManipInputPack& input_data, const std::vector<std::string> command)
namespace boxtools::actions{
  static std::unordered_map<std::string, void (*)(GroManipData&, const std::vector<std::string>&)> action_map;
  //use to initialize registry
  void merge(GroManipData&, const std::vector<std::string>&);
  void eraseres(GroManipData&, const std::vector<std::string>&);
  void keepres(GroManipData&, const std::vector<std::string>&);
  void duplicate(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyresname(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyresnameinv(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyatomname(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyatomnameinv(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyatomnameinvperiodic(GroManipData&, const std::vector<std::string>&);
  void trimvolumebyresnameinvperiodic(GroManipData&, const std::vector<std::string>&);
  void rotate(GroManipData&, const std::vector<std::string>&);
  void invrotate(GroManipData&, const std::vector<std::string>&);
  void rotate_vector(GroManipData&, const std::vector<std::string>&);
  void flip(GroManipData&, const std::vector<std::string>&);
  void translate(GroManipData&, const std::vector<std::string>&);
  void scale(GroManipData&, const std::vector<std::string>&);
  void loadgro(GroManipData&, const std::vector<std::string>&);
  void writegro(GroManipData&, const std::vector<std::string>&);
  void createprimitivevolume(GroManipData&, const std::vector<std::string>&);
  void createunionvolume(GroManipData&, const std::vector<std::string>&);
  void shrinkwrap(GroManipData&, const std::vector<std::string>&);
  void relabelatom(GroManipData&, const std::vector<std::string>&);
  void relabelres(GroManipData&, const std::vector<std::string>&);
  void respattern(GroManipData&, const std::vector<std::string>&);
  void setboxsize(GroManipData&, const std::vector<std::string>&);
  void outputmolecule(GroManipData&, const std::vector<std::string>&);
  void deleteoverlapping(GroManipData&, const std::vector<std::string>&);
  void pbccorrect(GroManipData&, const std::vector<std::string>&);
  void deleterandom(GroManipData&, const std::vector<std::string>&);
  //crystals
  void supercell(GroManipData&, const std::vector<std::string>&);
  void decoratefeldspar(GroManipData& data, const std::vector<std::string>& args);
  void hydrogenatefeldspar(GroManipData& data, const std::vector<std::string>& args);
  //use a mesh to trim the atoms
  void trimbymesh(GroManipData& data, const std::vector<std::string>& args);
  //for adding a single atom to the box, charge equalization
  void addatom(GroManipData& data, const std::vector<std::string>& args);
  //for getting atoms for union of spheres volume
  void printindicesnear(GroManipData& data, const std::vector<std::string>& args);
  //remove atoms with a certain number of nearest neighbors
  void hollow(GroManipData& data, const std::vector<std::string>& args);
  //wrap atoms to pbc
  void wrap(GroManipData& data, const std::vector<std::string>& args);
  //void outputbonds(GroManipData&, const std::vector<std::string>&);
  static inline void registerActions(){
    action_map.emplace("merge", &merge);
    action_map.emplace("eraseres", &eraseres);
    action_map.emplace("keepres", &keepres);
    action_map.emplace("duplicate", &duplicate);
    action_map.emplace("trimbyresname", &trimvolumebyresname);
    action_map.emplace("trimbyresnameinv", &trimvolumebyresnameinv);
    action_map.emplace("trimbyresnameinvp", &trimvolumebyresnameinvperiodic);
    action_map.emplace("trimbyatomname", &trimvolumebyatomname);
    action_map.emplace("trimbyatomnameinv", &trimvolumebyatomnameinv);
    action_map.emplace("trimbyatomnameinvp", &trimvolumebyatomnameinvperiodic);
    action_map.emplace("trimbymesh", &trimbymesh);
    action_map.emplace("rotate", &rotate);
    action_map.emplace("invrotate", &invrotate);
    action_map.emplace("rotate_vector", &rotate_vector);
    action_map.emplace("flip", &flip);
    action_map.emplace("translate", &translate);
    action_map.emplace("scale", &scale);
    action_map.emplace("loadgro", &loadgro);
    action_map.emplace("writegro", &writegro);
    action_map.emplace("volume", &createprimitivevolume);
    action_map.emplace("volumeunion", &createunionvolume);
    action_map.emplace("shrinkwrap", &shrinkwrap);
    action_map.emplace("relabelatom", &relabelatom);
    action_map.emplace("relabelres", &relabelres);
    action_map.emplace("setboxsize", &setboxsize);
    action_map.emplace("outputmolecule", &outputmolecule);
    action_map.emplace("respattern", &respattern);
    action_map.emplace("deleteoverlapping", &deleteoverlapping);
    action_map.emplace("pbccorrect", &pbccorrect);
    action_map.emplace("supercell", &supercell);
    action_map.emplace("decoratefeldspar", &decoratefeldspar);
    action_map.emplace("hydrogenatefeldspar", &hydrogenatefeldspar);
    action_map.emplace("printindicesnear", &printindicesnear);
    action_map.emplace("addatom", &addatom);
    action_map.emplace("deleterandom", &deleterandom);
    action_map.emplace("hollow", &hollow);
    action_map.emplace("wrap", &wrap);
    return;
  }
  static inline void performAction(GroManipData& data, std::vector<std::string> args){
    FANCY_ASSERT(args.size() >= 1, "Empty command issued to performAction()");
    std::string type = args[0];
    std::cout << type << std::endl;
    args.erase(args.begin());
    auto entry = action_map.find(type);
    FANCY_ASSERT(entry != action_map.end(), "Unable to find matching command");
    auto fptr = entry->second;
    (*fptr)(data, args);
    return;
  }
}