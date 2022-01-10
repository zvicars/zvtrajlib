#pragma once
#include "../datatypes.hpp"
namespace parse{
  //simple parser for ndx files, takes in a filename, finds and reads the file, constructs all atom groups
  IndexInfo parseNDX(std::string filename);
}