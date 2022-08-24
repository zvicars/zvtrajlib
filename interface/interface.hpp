#pragma once
#include <fstream>
#include <iostream>
#include "datatypes.hpp"
namespace xdr{ //to prevent any potential name clashes with the matrix type
#include "libxdr/xdrfile.h"
#include "libxdr/xdrfile_xtc.h"
#include "libxdr/xdrfile_trr.h"
}


class Trajectory{
public:
  virtual void getFrame(Box& box) = 0;
  virtual int nextFrame() = 0;

protected:
  int nframes_, frame_, natoms_, state_;
  float time_;
};


class XDRTrajectory : public Trajectory{
public:
  virtual void getFrame(Box& box) = 0;
  virtual int nextFrame() = 0;
protected:
  std::string trajfile_; //relevant files, only trajfile_ is strictly required
  xdr::XDRFILE * xd_;
  xdr::rvec * x_; //particle positions
  xdr::matrix box_; //box matrix
  float prec_;
};

class XTCTrajectory : public XDRTrajectory{
public:
  XTCTrajectory();
  XTCTrajectory(std::string trajfile);
  //frees positions array and closes file
  ~XTCTrajectory();
  //iterates to next frame in trajectory, sets all member variable info, returns 1 if no error, returns 0 if error, this value is also set in state_
  virtual int nextFrame();
  //replaces position and box vector data to correspond to current frame
  virtual void getFrame(Box& box);
};

class TRRTrajectory : public XDRTrajectory{
public:
  TRRTrajectory();
  TRRTrajectory(std::string trajfile);
  //frees positions array and closes file
  ~TRRTrajectory();
  //iterates to next frame in trajectory, sets all member variable info, returns 1 if no error, returns 0 if error, this value is also set in state_
  virtual int nextFrame();
  //replaces position and box vector data to correspond to current frame
  virtual void getFrame(Box& box);

protected:
  xdr::rvec * v_, * f_; //particle velocities
  float lambda_;
  int hasProp_;
};


void readNDX(std::string filename, Box& box);
void readGRO(std::string filename, Box& box);
void writeGRO_ov(std::string ofilename, const Box* box, std::array<double, 3> box_size, 
                 std::vector<int> indices, std::vector<std::array<double, 3> > positions);
void writeGRO(std::string ofilename, const Box* box);                
void writeXYZ_ov(std::string ofilename, const Box* box, std::array<double, 3> box_size, 
                 std::vector<int> indices, std::vector<std::array<double, 3> > positions);
void writeXYZ_ov_append(std::string ofilename, std::array<double, 3> box_size, 
                 std::vector<std::string> names, std::vector<std::array<double, 3> > positions);
inline void readTOP(std::string filename, Box& box){ return; }


class GROTrajectory : public Trajectory{
public:
  GROTrajectory(std::string grofile){
    grofile_ = grofile;
    return;
  }
  virtual void getFrame(Box& box){
    box.time = 0.0;
    box.frame = 1;
    readGRO(grofile_, box);
  }
  int nextFrame(){
    if(firstFrame){
      firstFrame = 0;
      return 1;
    }
    return 0;
  };
private:
  bool firstFrame = 1;
  std::string grofile_;
};

static inline Trajectory* loadTrajectory(std::string trajectory_file){
  std::string filetype = trajectory_file.substr(trajectory_file.rfind('.')+1);
  if(filetype == "xtc"){
    return new XTCTrajectory(trajectory_file);
  }
  if (filetype == "trr"){
    std::cout << "trr files have not been implemented yet!" << std::endl;
    throw 1;
    return 0;
  }
  if (filetype == "gro"){
    return new GROTrajectory(trajectory_file);
  } 
  std::cout << "Unrecognized trajectory filetype, " << filetype << ".";
  throw 1;
  return 0;
}