#include "parseNDX.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include "string_ops.hpp"
IndexInfo parseNDX(std::string filename){

  std::map<std::string, Vec<int> > output_map;

  std::ifstream ifile(filename);
  if(!ifile.is_open()){
    std::cout << "Failed to open .ndx file" << std::endl;
    throw 0; 
  }
  int group_counter = 0;
  std::vector<std::string> groups;
  std::vector<int> index_counts;
  std::string line;
  while(std::getline(ifile, line)){
    if(line.length() == 0) continue;
    if(line.at(0) == ';') continue;
    int pos1 = line.find("[");
    int pos2 = line.find("]");
    if(pos1 == std::string::npos) continue;
    if(pos2 == std::string::npos || pos2 < pos1) continue;
    //if I get this far I found a valid index label
    std::string label = line.substr(pos1+1, pos2-pos1-1);
    trim(label);
    std::vector<int> indexes;
    //read the label, should be nothing but an array of numbers until another line with brackets
    auto lastpos = ifile.tellg();
    int counter = 0;
    while(std::getline(ifile, line)){  
      counter++;
      if(line.size() == 0) continue;
      if(line.at(0) == ';') continue;
      //leave the loop once a line with a square bracket is found
      int pos1 = line.find("[");
      if(pos1 != std::string::npos){
        ifile.seekg(lastpos);
        break;
      }
      std::stringstream ss(line);
      int number;
      while(ss >> number){
        indexes.push_back(number);
      }
      //store previous line's location so the outer while loop gets the correct line once this loop exits
      lastpos = ifile.tellg();
    }
    output_map.insert(std::pair<std::string, Vec<int> >(label, indexes));
    groups.push_back(label);
    index_counts.push_back(indexes.size());
    group_counter++;
  }
  ifile.close();

  std::cout << "Found " << group_counter << " groups. In index file " << filename << "\n";
  for(int i = 0; i < groups.size(); i++){
    std::cout << groups[i] << "  " << index_counts[i] << " atoms." << "\n";
  }
  std::cout << std::endl;
  IndexInfo a1;
  a1.indexes = output_map;
  return a1;
}