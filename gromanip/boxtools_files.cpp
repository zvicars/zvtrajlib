#include "boxtools.hpp"
#include <algorithm>
#include "../tools/StringTools.hpp"
namespace boxtools{
  std::string loadFileIntoString(std::string filename){
    std::ifstream ifile(filename);
    FANCY_ASSERT(ifile.is_open(), "Failed to open input file in loadFileIntoString()");
    std::string line;
    std::stringstream buffer;
    while(std::getline(ifile, line)){
      buffer << line << "\n";        
    }
    return buffer.str();
  }
  std::string getGMXTaggedParams(std::string filedata, std::string key){
    std::stringstream ifile(filedata);
    std::string line;
    std::string retval = "";
    bool keyFound = 0; //has the key been found/is it reading data currently
    while(std::getline(ifile, line)){
      line = StringTools::trimWhitespace(line);
      if(line.size() == 0) continue;
      if(line.at(0) == ';' || line.at(0) == '#') continue;
      if(keyFound == 1){
        if(line.size() == 0) break;
        if(line.at(0) == '[') break;
        retval = retval + line + "\n";
      }
      if(line == "[ " + key + " ]"){
        keyFound = 1;
      }
    }
    return retval;
  }
}