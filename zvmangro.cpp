#include "gromanip/actions.hpp"
#include <fstream>
#include <string>
#include <sstream>
int main(int argc, char **argv)
{
  boxtools::actions::registerActions();
  FANCY_ASSERT(argc == 2, "Analysis code only accepts a single input that specifies the op input file.");
  std::string input_file = argv[1];
  std::ifstream ifile(input_file);
  FANCY_ASSERT(ifile.is_open(), "Failed to open input file stream for zvmangro exectuable");
  GroManipData data;
  std::string line;
  while(std::getline(ifile, line)){
    if(line.at(0) == '#' || line.length() == 0 || line.at(0) == ';' || line.at(0) == ' ' || line.length()==0) continue;
    std::stringstream ss(line);
    std::vector<std::string> line_args;
    std::string lineval;
    while(ss >> lineval){
      line_args.push_back(lineval);
    }
    boxtools::actions::performAction(data, line_args);
  }
  std::cout << "done with all commands..." << std::endl;
  return 0;
}