#include "Calc_Write_Xtc.hpp"

Calc_Write_Xtc::Calc_Write_Xtc(InputPack& input) : Calculation{input} {
  filename_ = name_ + ".xtc";
  input.params().readString("filename", KeyType::Optional, filename_);

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  initialized_ = 0;
  return;
}

void Calc_Write_Xtc::update(){
  Calculation::update();

  xdr_time_ = current_time_;
  xdr_step_ = current_frame_;
  xdr_natoms_ = box->atoms.size();

  if(!initialized_){
    xdr_prec_ = 1000;
    xdr_x_ = new xdr::rvec[xdr_natoms_];
    output_handle_ = xdr::xdrfile_open(filename_.c_str(), "w");
    FANCY_ASSERT(output_handle_!=0, "xtc output file is null");
    initialized_ = 1;
  }

  for(int i = 0; i < box->atoms.size(); i++){
    for(int j = 0; j < 3; j++){
      xdr_x_[i][j] = -1.0f;
    }
  }

  auto indices = atom_group_->getIndices();
  for(int i = 0; i < indices.size(); i++){
    int index = indices[i];
    for(int j = 0; j < 3; j++){
      xdr_x_[index][j] = box->atoms[index].x[j];
    }
  }

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
        xdr_box_[i][j] = box->boxvec[i][j];
    }
  }

  return;
}

void Calc_Write_Xtc::output(){
  if(!doOutput()) return;
  xdr::write_xtc(output_handle_, xdr_natoms_, xdr_step_, xdr_time_, xdr_box_, xdr_x_, xdr_prec_);
  return;
}

void Calc_Write_Xtc::finalOutput(){
  //maybe this would be better in the destructor?
  xdr::xdrfile_close(output_handle_);
  delete xdr_x_;
  return;
}