#include "gromanip/actions.hpp"

int main(int argc, char **argv)
{
  FANCY_ASSERT(argc == 2, "Analysis code only accepts a single input that specifies the op input file.");
  std::string input_file = argv[1];
  InputParser input_parser;
  //load json with key/value pairs
  ParameterPack master_pack = input_parser.parseFile(input_file); 
  using KeyType = ParameterPack::KeyType;

  //output gro file
  std::string output_file;
  master_pack.readString("output_file", KeyType::Required, output_file);

  //


}