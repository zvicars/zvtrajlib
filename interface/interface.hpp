#pragma once
#include <fstream>
#include "datatypes.hpp"
#include "boxtools.hpp"
namespace xdr{ //to prevent any potential name clashes with the matrix type
#include "libxdr/xdrfile.h"
#include "libxdr/xdrfile_xtc.h"
}
class XDRTrajectory{
public:
  XDRTrajectory();
  XDRTrajectory(std::string trajfile);
  //frees positions array and closes file
  ~XDRTrajectory();
  //iterates to next frame in trajectory, sets all member variable info, returns 1 if no error, returns 0 if error, this value is also set in state_
  int nextFrame();
  //replaces position and box vector data to correspond to current frame
  void getFrame(Box& box);

private:
  std::string trajfile_; //relevant files, only trajfile_ is strictly required
  xdr::XDRFILE * xd_;
  xdr::rvec * x_; //particle positions
  xdr::matrix box_; //box matrix
  int nframes_, frame_, natoms_, state_;
  float time_, prec_;
};

void readNDX(std::string filename, Box& box);
void readGRO(std::string filename, Box& box);
void writeGRO_ov(std::string ofilename, const Box* box, std::array<double, 3> box_size, 
                 std::vector<int> indices, std::vector<std::array<double, 3> > positions);
void writeGRO(std::string ofilename, const Box* box);                
void writeXYZ_ov(std::string ofilename, const Box* box, std::array<double, 3> box_size, 
                 std::vector<int> indices, std::vector<std::array<double, 3> > positions);
inline void readTOP(std::string filename, Box& box){
  return;
}
