// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "StringTools.hpp"
namespace StringTools {

bool stringToBool(const std::string& str)
{
  std::string lowercase_str = toLowercase(str);

  if ( lowercase_str == "yes" || lowercase_str == "true" || 
       lowercase_str == "1"   || lowercase_str[0] == 'y' ) { 
    return true; 
  }
  else if ( lowercase_str == "no" || lowercase_str == "false" ||
            lowercase_str == "0"  || lowercase_str[0] == 'n' ) { 
    return false;
  }
  else {
    FANCY_ASSERT( false,
      "unable to interpret '" << str << "' as as true/yes/1 or false/no/0: check your input." );
  }
}


std::string toLowercase(const std::string& str)
{
  std::string lowercase_str = str;
  std::transform(lowercase_str.begin(), lowercase_str.end(), lowercase_str.begin(), ::tolower);
  return lowercase_str;
}


std::string getFileExtension(const std::string& file)
{
  auto last_period = file.find_last_of(".");
  if ( last_period != std::string::npos ) {
    return file.substr(last_period + 1);
  }
  else {
    return "";
  }
}


std::string trimWhitespace(const std::string& str, const std::string& whitespace)
{
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos) { return ""; }

  const auto strEnd   = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}


void removeTrailingCommentFromLine(
  const std::string& str, const std::string& comment_chars,
  std::string& trimmed_str, std::string& comment
)
{
  auto pos_comment = str.find_first_of(comment_chars);

  if ( pos_comment != std::string::npos ) {
    if ( pos_comment > 0 ) {
      auto trimmed_str_length = pos_comment - 1;
      trimmed_str = str.substr(0, trimmed_str_length);
    }
    else {
      // Comment is the entire line
      trimmed_str = "";
    }
    comment = str.substr(pos_comment + comment_chars.size());
  }
  else {
    // No comment
    trimmed_str = str;
    comment     = "";
  }
}


bool isNumber(const std::string& str)
{
  std::stringstream ss( str );
  float value;
  ss >> value;
  return ( ! ss.fail() );
}


// Template specializations
template<>
std::string stringToValue<std::string>(const std::string& str) {
  return str;
}

template<>
bool stringToValue<bool>(const std::string& str) {
  return stringToBool(str);
}

void swap(std::string& s1, std::string& s2){
  std::string temp = s1;
  s1 = s2;
  s2 = temp;
  return;
}

}  // end namespace StringTools