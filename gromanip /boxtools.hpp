//takes a box and ensures that resnumbers and atomnumbers range from 1->NRES and 1->NATOMS respectively
#pragma once
#include "../interface/datatypes.hpp"
#include "volume.hpp"
namespace boxtools{
  bool renumberBox(Box& box);
  bool deleteRestype(Box& box, std::string resname);
  Box mergeBox(const Box& b1, const Box& b2);
  void removeResNumbers(Box& box, std::vector<int> res_list);
  std::vector<int> getResnrWithinVolume(const Box& box, const Volume& volume);
  std::vector<int> getResnrWithinVolumebyAtomName(const Box& box, const Volume& volume, std::string at_name);
  std::vector<int> getResnrAllButResname(const Box& box, const Volume& volume, std::string res_name);
  void rotateEulerAngles(Box& box, Vec3<double> angles);
  void translateAtoms(Box& box, Vec3<double> offset);
  void shrinkWrap(Box& box);
  bool checkBoxValid(const Box& box){
    if(box.atoms.size() > 0 && box.hasNamedAtoms) return 1;
    else throw 1;
    return 0; 
  }
}