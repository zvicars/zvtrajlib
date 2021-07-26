#include "parseGRO.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include "string_ops.hpp"

//1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
//0-4 RESNR
//5-9 RESNAME
//10-14 ATOMNAME
//15-19 ATOMNUMBER
//20-27 POSX
//28-35 POSY
//36-43 POSZ
//44-51 VX
//52-59 VY
//60-67 VZ
void parseGRO(std::string filename, Box& box){
  std::ifstream ifile(filename);
  if(!ifile.is_open()){
    std::cout << "Failed to open ifile" << std::endl;
    throw 0;
  }

  std::string line;
  std::getline(ifile, line); //comment line
  std::getline(ifile, line);
  int natoms = std::stod(line); //number of atoms
  box.atoms.resize(natoms);
  for(int i = 0; i < natoms; i++ ){
    std::getline(ifile, line); //these should all be normal grofile lines in fixed-column format
    int resnr = std::stoi(line.substr(0, 5));
    std::string resname = line.substr(5, 5);
    trim(resname);
    std::string atomname = line.substr(10, 5);
    int atomnumber = std::stoi(line.substr(15, 5));
    Vec3<double> position, velocity;
    position[0] = std::stod(line.substr(20, 8));
    position[1] = std::stod(line.substr(28, 8));
    position[2] = std::stod(line.substr(36, 8));
    if(line.length() >= 68){
      velocity[0] = std::stod(line.substr(44, 8));
      velocity[1] = std::stod(line.substr(52, 8));
      velocity[2] = std::stod(line.substr(60, 8));
    }
    else{
      velocity[0] = 0.0;
      velocity[1] = 0.0;
      velocity[2] = 0.0;
    }
    box.atoms[i].name = atomname;
    box.atoms[i].x = position;
    box.atoms[i].v = velocity;
    box.atoms[i].resname = resname;
    box.atoms[i].resnr = resnr;
    box.atoms[i].index = atomnumber;
  }
  box.time = 0;
  box.frame = 0;
  box.frame_counter++;
  box.hasNamedAtoms = 1;
  return;
}