#include "interface/interface.hpp"
#include "io/InputParser.hpp"
#include "analysis/parser.hpp" 
#include <iostream>
int main()
{

  InputParser input_parser;
  auto input_pack = input_parser.parseFile(op_input_file_);




  op_input_file_ = "test.input";
  InputParser input_parser;
  using KeyType = ParameterPack::KeyType;

  auto input_pack = input_parser.parseFile(op_input_file_);
  std::string trajectory_file, index_file, topology_file;
  input_pack.readString("trajectory", KeyType::Required, trajectory_file_);
  auto found_ndx = input_pack.readString("index", KeyType::Optional, index_file);
  auto found_top = input_pack.readString("topology", KeyType::Optional, topology_file);

  std::cout << "Starting test..." << std::endl;
  XDRTrajectory traj(trajectory_file);
  if(found_ndx) traj.addIndexFile(index_file);
  if(found_top) traj.addTopFile(top_file);
  Box b1;



  while(traj.nextFrame()){
    traj.getFrame(b1);
    double pvparms[6] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    int count = OP_Nv(pvparms, b1);
    std::cout << b1.time << "  " << count << std::endl;
  }
  return 0;
}