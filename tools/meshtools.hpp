#pragma once
#include "meshtools_files/meshfiles.hpp"
#include "meshtools_files/obj.hpp"
#include "StringTools.hpp"
static inline Mesh getMesh(std::string filename){
  
  Mesh meshout;

  std::ifstream t(filename);
  if(!t.is_open()){
    std::cout << "Failed to open file." << std::endl;
    throw 1;
  }
  std::string str((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());
  t.close();

  std::string file_ending = filename.substr(filename.rfind('.') + 1);
  file_ending = StringTools::trimWhitespace(file_ending);
  if(file_ending == "obj"){
    ObjFile obj;
    obj.load(str);
    meshout = obj.getMesh();
  }
  else{
    std::cout << "unrecognized filetype" << std::endl;
    throw 2;
  }

  return meshout;
}

bool MollerTrombore(std::array<double,3> rayOrigin, 
                           std::array<double,3> rayVector, 
                           std::array<std::array<double,3>, 3> inTriangle,
                           std::array<double,3>& outIntersectionPoint);