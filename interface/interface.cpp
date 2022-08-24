#include "interface.hpp"
#include "fileparsers/parseNDX.hpp"
#include "fileparsers/parseGRO.hpp"
#include <iostream>
#include <cstring>

void readNDX(std::string filename, Box& box){
  box.idxinfo = parse::parseNDX(filename);
  box.hasIndexes = 1;
  return;
}
void readGRO(std::string filename, Box& box){
  parse::parseGRO(filename, box);
  return;
}
void writeGRO_ov(std::string ofilename, const Box* box, std::array<double, 3> box_size, 
                 std::vector<int> indices, std::vector<std::array<double, 3> > positions)
{
  parse::writeGRO_override(ofilename, box, box_size, indices, positions);
  return;
}
void writeGRO(std::string ofilename, const Box* box)
{
  parse::writeGRO(ofilename, box);
  return;
}

void writeXYZ_ov(std::string ofilename, const Box* box, std::array<double, 3> box_size, 
                 std::vector<int> indices, std::vector<std::array<double, 3> > positions)
{
  std::ofstream ofile(ofilename);
  if(!ofile.is_open()){
    std::cerr << "Failed to open output file stream when writing xyz file." << std::endl;
    throw;
  }

  ofile << indices.size() << "\n\n";
  for(int i = 0; i < indices.size(); i++){
    auto pos = positions[i];
    auto name = box->atoms[indices[i]].name;
    ofile << name << "     " << 10.0*pos[0] << "     " << 10.0*pos[1] << "     " << 10.0*pos[2] << "\n";
  }
  ofile << std::endl;
  return;
}

void writeXYZ_ov_append(std::string ofilename, std::array<double, 3> box_size, 
                 std::vector<std::string> names, std::vector<std::array<double, 3> > positions)
{
  std::ofstream ofile(ofilename, std::ios_base::app);
  if(!ofile.is_open()){
    std::cerr << "Failed to open output file stream when writing xyz file." << std::endl;
    throw;
  }
  ofile << names.size() << "\n\n";
  for(int i = 0; i < names.size(); i++){
    auto pos = positions[i];
    auto name = names[i];
    ofile << name << "     " << 10.0*pos[0] << "     " << 10.0*pos[1] << "     " << 10.0*pos[2] << "\n";
  }
  ofile << std::endl;
  return;
}