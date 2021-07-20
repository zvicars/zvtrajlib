#include "interface.hpp"
#include <iostream>
#include <cstring>
XDRTrajectory::XDRTrajectory(){
  //default constructor will set some invalid values that must be set later
  xdr::XDRFILE * xd = 0;
  xdr::rvec * x = 0; //particle positions
  nframes_ = -1;
  return;
}

XDRTrajectory::XDRTrajectory(std::string trajfile){
  trajfile_= trajfile;
  char tf[1024];
  strcpy(tf, trajfile.c_str());
  int result = xdr::read_xtc_natoms(tf, &natoms_);
  if(result == xdr::exdrFILENOTFOUND){
    std::cout << "File not found." << std::endl;
    throw result;
    return;
  }
	x_ = new xdr::rvec[natoms_];
	xd_ = xdr::xdrfile_open(tf, "r");
  nframes_ = 0;
  return;
}

XDRTrajectory::~XDRTrajectory(){
  delete[] x_;
  xdr::xdrfile_close(xd_);
}

int XDRTrajectory::nextFrame(){
  state_ = xdr::read_xtc(xd_, natoms_, &frame_, &time_, box_, x_, &prec_);
  if(state_ != xdr::exdrOK){ // general error if frame fails to advance, to prevent having to care about xdrfile stuff outside of class
    return 0;
  }
  nframes_++;
  return 1; //successful frame read
}

void XDRTrajectory::getFrame(Box& box){
  box.atoms.resize(natoms_);
  for(int i = 0; i < natoms_; i++)
  {
    box.atoms[i].x[0] = x_[i][0];
    box.atoms[i].x[1] = x_[i][1];
    box.atoms[i].x[2] = x_[i][2];
  }
  for(int i = 0; i < 3; i++)
  for(int j = 0; j < 3; j++)
  {
    box.boxvec[i][j] = box_[i][j];
  }
  box.time = time_;
  box.frame = frame_;
  box.frame_counter = nframes_-1;
  return;
}