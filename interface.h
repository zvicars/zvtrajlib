#pragma once
#include <string>
#include <cstring>
#include <fstream>
#include "datatypes.h"
namespace xdr{ //to prevent any potential name clashes with the matrix type
#include "libxdr/xdrfile.h"
#include "libxdr/xdrfile_xtc.h"
}
class XDRTrajectory{
public:
  XDRTrajectory();
  XDRTrajectory(std::string trajfile);
  ~XDRTrajectory();
  int nextFrame();
  void getFrame(Box& box);

private:
  std::string trajfile_; //relevant files, only trajfile_ is strictly required
  xdr::XDRFILE * xd_;
  xdr::rvec * x_; //particle positions
  xdr::matrix box_; //box matrix
  int nframes_, frame_, natoms_, state_;
  float time_, prec_;
};