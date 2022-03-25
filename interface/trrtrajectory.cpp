#include "interface.hpp"
#include "fileparsers/parseNDX.hpp"
#include "fileparsers/parseGRO.hpp"
#include <iostream>
#include <cstring>
TRRTrajectory::TRRTrajectory(){
  //default constructor will set some invalid values that must be set later
  xdr::XDRFILE * xd = 0;
  xdr::rvec * x = 0; //particle positions
  nframes_ = -1;
  return;
}

TRRTrajectory::TRRTrajectory(std::string trajfile){
  trajfile_= trajfile;
  char tf[1024];
  strcpy(tf, trajfile.c_str());
  int result = xdr::read_trr_natoms(tf, &natoms_);
  if(result == xdr::exdrFILENOTFOUND){
    std::cout << "File not found." << std::endl;
    throw result;
    return;
  }
	x_ = new xdr::rvec[natoms_];
  v_ = new xdr::rvec[natoms_];
	xd_ = xdr::xdrfile_open(tf, "r");
  nframes_ = 0;
  return;
}

TRRTrajectory::~TRRTrajectory(){
  delete[] x_;
  xdr::xdrfile_close(xd_);
}

int TRRTrajectory::nextFrame(){
  state_ = xdr::read_trr(xd_, natoms_, &frame_, &time_, &lambda_, box_, x_, v_, f_, &hasProp_);
  if(state_ != xdr::exdrOK){ // general error if frame fails to advance, to prevent having to care about xdrfile stuff outside of class
    return 0;
  }
  nframes_++;
  return 1; //successful frame read
}

void TRRTrajectory::getFrame(Box& box){
  if(box.atoms.size() == 0) box.atoms.resize(natoms_); //if data hasn't been loaded in already, ensure that there's enough space allocated
  if(box.atoms.size() != natoms_){
    std::cout << "Mismatch between box atoms and TRR atoms, make sure all of your trajectory files are coming from the right simulation!" << std::endl;
    throw 0;
  }
  for(int i = 0; i < natoms_; i++)
  {
    box.atoms[i].x[0] = x_[i][0];
    box.atoms[i].x[1] = x_[i][1];
    box.atoms[i].x[2] = x_[i][2];
    box.atoms[i].v[0] = v_[i][0];
    box.atoms[i].v[1] = v_[i][1];
    box.atoms[i].v[2] = v_[i][2];    
    box.atoms[i].f[0] = f_[i][0];
    box.atoms[i].f[1] = f_[i][1];
    box.atoms[i].f[2] = f_[i][2];        
  }
  for(int i = 0; i < 3; i++)
  for(int j = 0; j < 3; j++)
  {
    box.boxvec[i][j] = box_[i][j];
  }
  box.time = time_;
  box.frame = frame_;
  box.frame_counter++;
  box.hasForces = 1;
  box.hasVelocities = 1;
  return;
}