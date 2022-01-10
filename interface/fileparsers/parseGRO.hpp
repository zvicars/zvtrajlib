#include "../datatypes.hpp"
namespace parse{
  //gro files contain a single frame and atom names, will load that single frame
  void parseGRO(std::string filename, Box& box);
  //write a gro file using information from the box, but replacing the positions and box_vec with manually specified values
  void writeGRO_override(std::string ofilename, const Box* box, std::array<double, 3> box_size, 
                        std::vector<int> indices, std::vector<std::array<double, 3> > positions);
  void writeGRO(std::string ofilename, const Box* box);
}