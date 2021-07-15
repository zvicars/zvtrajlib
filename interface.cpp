#include "interface.h"
XDRTrajectory::XDRTrajectory(){
  //default constructor will set some invalid values that must be set later
  xdr::XDRFILE * xd = 0;
  xdr::rvec * x = 0; //particle positions
  nframes_ = -1;
  return;
}

XDRTrajectory::XDRTrajectory(std::string trajfile){
  trajfile_= trajfile;
  char * tf = 0;
  tf = strdup(trajfile.c_str());
  xdr::read_xtc_natoms(tf, &natoms_);
	x_ = new xdr::rvec[natoms_];
	xd_ = xdr::xdrfile_open(tf, "r");
  nframes_ = 0;
  free(tf);
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
  return;
}