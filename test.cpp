#include "interface/interface.hpp"
#include "analysis/ProbeVolume.hpp"
#include "analysis/Calculation.hpp"
int main()
{
  std::string op_input_file_ = "test.input";
  std::string trajectory_file_;
  InputParser input_parser;
  auto input_pack = input_parser.parseFile(op_input_file_);
  using KeyType = ParameterPack::KeyType;
  input_pack.readString("trajectory", KeyType::Required, trajectory_file_);

  std::cout << "Starting test..." << std::endl;
  XDRTrajectory traj(trajectory_file_);
  Box b1;
  while(traj.nextFrame()){
    traj.getFrame(b1);
  }
  return 0;
}