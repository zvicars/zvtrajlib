//takes a box and ensures that resnumbers and atomnumbers range from 1->NRES and 1->NATOMS respectively
#pragma once
#include "../interface/datatypes.hpp"
#include "boxtools_datatypes.hpp"
#include "volume.hpp"

namespace boxtools{
  template <class T> 
  struct ParameterTable{
    public:
      T& operator[](int i){
        return entries[i];
      }
      void push_back(const T& input){
        entries.push_back(input);
      }
      int size(){
        return entries.size();
      }
      std::vector<T> entries;
  };

  template <class T>
  static ParameterTable<T> generateTable(std::string input_data){
    ParameterTable<T> output;
    std::stringstream ss(input_data);
    std::string line;
    while(std::getline(ss,line)){
      T instance;
      instance.read(line);
      output.push_back(instance);
    }
    return output; 
  }


  std::string loadFileIntoString(std::string filename);
  std::string getGMXTaggedParams(std::string filedata, std::string key);
  std::vector<AtomInst> makeAtoms(const Box& box, const ParameterTable<BTAtomType>& atable);
  std::vector<PosResInst> makeRestraints(const Box& box, const ParameterTable<PosResType>& t);
  std::vector<BondInst> makeBonds(const Box& box, const ParameterTable<BondType>& btable);
  std::vector<AngleInst> makeAngles(const Box& box, const ParameterTable<AngleType>& atable);
  std::vector<AngleInst> makeAnglesUsingBonds(const Box& box, const ParameterTable<AngleType>& atable, const std::vector<BondInst>& bonds);
  std::vector<BondInst> makePeriodicBonds(Box& b1, const ParameterTable<PBCBondType>& t);
  std::vector<ExclusionInst> makeExclusions(Box& b1, const ParameterTable<ExclusionType>& t);
  std::vector<ConstraintInst> makeConstraints(Box& b1, const ParameterTable<ConstraintType>& t);
  std::vector<Vsite3Inst>  makeVsites(Box& b1, const ParameterTable<Vsite3Type>& t);

  bool renumberBox(Box& box);
  bool renumberBoxSeq(Box& box);
  bool deleteRestype(Box& box, std::string resname);
  bool deleteNotRestype(Box& box, std::string resname);
  Box mergeBox(const Box& b1, const Box& b2);
  void removeResNumbers(Box& box, std::vector<int> res_list);
  std::vector<int> getResnrWithinVolume(const Box& box, const Volume& volume);
  std::vector<int> getResnrNotWithinVolume(const Box& box, const Volume& volume);
  std::vector<int> getResnrWithinVolumebyAtomName(const Box& box, const Volume& volume, std::string at_name);
  std::vector<int> getResnrNotWithinVolumebyAtomName(const Box& box, const Volume& volume, std::string at_name);
  std::vector<int> getResnrAllButResname(const Box& box, const Volume& volume, std::string res_name);
  void rotateEulerAngles(Box& box, Vec3<double> angles);
  void invrotateEulerAngles(Box& box, Vec3<double> angles);
  void translateAtoms(Box& box, Vec3<double> offset);
  void scaleAtoms(Box& box, Vec3<double> offset);
  void wrapPBC(Box& box);
  void shrinkWrap(Box& box, double buffer = 0.0);
  void relabelAtom(Box& box, std::string original_name, std::string new_name);
  void relabelRes(Box& box, std::string original_name, std::string new_name);
  void setBoxSize(Box& box, Vec3<double> newsize);
  void rotateVectorCOM(Box& box, Vec3<double> v1, Vec3<double> v2);
  bool setAtomParams(Box& box, const ParameterTable<AtomType>& atable);
  static inline bool checkBoxValid(const Box& box){
    if(box.atoms.size() > 0 && box.hasNamedAtoms) return 1;
    else throw 1;
    return 0; 
  }
}