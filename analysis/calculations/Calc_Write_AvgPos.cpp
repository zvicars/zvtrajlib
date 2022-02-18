#include "Calc_Write_AvgPos.hpp"
#include "../helper/pbc_functions.hpp"
Calc_Write_AvgPos::Calc_Write_AvgPos(InputPack& input) : Calculation{input} {
  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  input.params().readString("output_file", KeyType::Required, ofilename_);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  FANCY_ASSERT(!atom_group_->isDynamic(), "Calc_Write_AvgPos requires the atom group to be static.");
  initialized_ = 0;
  frame_counter_ = 0;
  avg_box_size_ = {0.0};
  return;
}

void Calc_Write_AvgPos::update(){
  if(hasUpdated()) return;
  Calculation::update();
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  if(!initialized_){
    natoms_ = atom_group_->getIndexCount();
    setRefPosArray(); //set initial atom positions to ensure correct pbc is used for edge atoms, also stores indices
    initialized_ = 1;
    FANCY_ASSERT(box->hasNamedAtoms, "Calc_Write_AvgPos needs named atoms for gro file output.");
  }
  return;
}
void Calc_Write_AvgPos::setRefPosArray(){

  RefIdx_.resize(natoms_);
  RefPos_.resize(natoms_);
  ComPos_.resize(natoms_, {0.0});
  for(int i = 0; i < natoms_; i++){
    int index = atom_group_->getIndex(i);
    RefIdx_[i] = index;
    RefPos_[i] = box->atoms[index].x;
  }
  return;
}

void Calc_Write_AvgPos::calculate(){
  if(!doCalculate()) return;
  #pragma omp parallel for
  for(int i = 0; i < natoms_; i++){
    auto newpos = getNearestPeriodicPosition(box->atoms[RefIdx_[i]].x, RefPos_[i], box_size_);
    for(int j = 0; j < 3; j++){
      ComPos_[i][j] += newpos[j];
    }
  }
  for(int j = 0; j < 3; j++){
    avg_box_size_[j] += box_size_[j];
  } 
  frame_counter_++;
  return;
}

void Calc_Write_AvgPos::output(){
  if(!doOutput()) return;
  return;
}

void Calc_Write_AvgPos::finalOutput(){
  double invcount = 1.0/(double)frame_counter_;
  for(int i = 0; i < natoms_; i++){
    for(int j = 0; j < 3; j++){
      ComPos_[i][j] *= invcount;
    }
  }
  for(int j = 0; j < 3; j++){
    avg_box_size_[j] *= invcount;
  }
  writeGRO_ov(ofilename_, box, avg_box_size_, RefIdx_, ComPos_);
  //writeXYZ_ov(ofilename_, box, avg_box_size_, RefIdx_, ComPos_);
  return;
}