// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef STRING_TOOLS_H
#define STRING_TOOLS_H

#include <algorithm>
#include <array>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Assert.hpp"

// Helper functions for dealing with std::strings
namespace StringTools {


//----------------------------//
//----- Sanitizing Input -----//
//----------------------------//

// Removes leading and trailing whitespace
// - Solution of StackOverflow user "GManNickG" (Nov 25 '09 at 16:29)
// - See: https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
std::string trimWhitespace(
  const std::string& str,
  const std::string& whitespace = " \t" // defines "whitespace"
);

// Trim the comment off the end of a string, where the beginning of the comment is
// *any* of the characters in "comment_chars" (e.g. the character #)
void removeTrailingCommentFromLine(
  const std::string& str,
  const std::string& comment_chars,
  std::string& trimmed_str,
  std::string& comment
);


//--------------------------//
//----- Tokenizization -----//
//--------------------------//

// Split whitespace-delimited string into tokens
template<typename T = std::string>
std::vector<T> split(const std::string& str);

// Split strings delimited by 'delim' into tokens
template<typename T = std::string>
std::vector<T> split(const std::string& str, const char delim);



//-----------------------//
//----- Conversions -----//
//-----------------------//

// Convert a std::string to a value of the given type
template<typename T> 
T stringToValue(const std::string& str);

// Converts "yes/no", "true/false", etc. to the appropriate bool
bool stringToBool(const std::string& str);

template<typename T, std::size_t dim>
std::array<T,dim> stringsToArray(const std::vector<std::string>& strings);

template<typename T>
std::vector<T> stringsToVector(const std::vector<std::string>& strings);



//------------------//
//----- Output -----//
//------------------//

// Helper class for converting a value to a fixed-width string
template<typename T>
class FixedWidthValue
{
 public:
  FixedWidthValue(const T& value, const int width = -1):
    value_(value), width_(-1)
  {}

  const T&  getValue() const { return value_; }
  const int getWidth() const { return width_; }

  private:
  const T& value_;
 int width_;
};

// Prints a fixed-width string to the indicated stream
// - Assumes any desired formatting flags have already been set (e.g. by std::fixed)
// - The string is truncated as necessary according to the stream's 'adjustfield' flag
//   - std::left     - trim from the right
//   - std::right    - trim from the left
//   - std::internal - unsupported (throws an exception)
//   - (other)       - same as std::right
// - If 'width' is non-negative, calls os.setw(width)
//   - Else uses whatever value (if any) is already present in the stream
template<typename T>
std::ostream& operator<<(std::ostream& os, const FixedWidthValue<T> value);



//----------------//
//----- Misc -----//
//----------------//

// Convert a string to all lowercase
std::string toLowercase(const std::string& str);

// Returns the file extension of the given file name (if it can be found),
// else it returns an empty string
// - File extension is everything from the last period (.) to the end
// - e.g. if file = "myFile.dat", returns "dat"
//        if file = "myFile", returns ""
// - TODO: move to FileSystem?
std::string getFileExtension(const std::string& file);

// Returns true if the string can be parsed as a number
bool isNumber(const std::string& str);





//---------------------------//
//----- Implementations -----//
//---------------------------//





template<typename T>
std::vector<T> split(const std::string& str)
{
  std::vector<T> tokens;
  std::stringstream ss(str);
  T token;
  while ( ss >> token ) {
    tokens.push_back(token);
  }

  return tokens;
}


template<typename T>
std::vector<T> split(const std::string& str, const char delim)
{
  std::vector<T> tokens;
  std::stringstream ss(str);

  std::string line;
  while ( getline(ss, line, delim) ) {
    tokens.push_back( stringToValue<T>(line) );
  }

  return tokens;
}

template<typename T> 
T stringToValue(const std::string& str) {
  std::stringstream ss( str );
  T value;
  ss >> value;

  FANCY_ASSERT( ! ss.fail(), "unable to convert '" << str << "' to the desired type" );

  return value;
}


template<typename T, std::size_t dim>
std::array<T,dim> stringsToArray(const std::vector<std::string>& strings) {
  FANCY_ASSERT( strings.size() == dim, "size mismatch" );

  std::array<T,dim> arr;
  std::transform( strings.begin(), strings.end(), arr.begin(),
    [](const std::string& s) { return StringTools::stringToValue<T>(s); } );
  return arr;
}

template<typename T>
std::vector<T> stringsToVector(const std::vector<std::string>& strings) {
  std::vector<T> vec( strings.size() );
  std::transform( strings.begin(), strings.end(), vec.begin(),
    [](const std::string& s) { return StringTools::stringToValue<T>(s); } );
  return vec;
}

//shorter comma-delimited variant, [x,y,z]
template<typename T>
std::vector<T> stringToVector(const std::string& str) {
  int pos1 = str.find('[');
  int pos2 = str.rfind(']');
  std::string str_data = str.substr(pos1+1, pos2-pos1-1);
  std::stringstream ss(str_data);
  std::vector<T> vec;
  std::string token;
  while(std::getline(ss, token, ',')) {
      std::stringstream ss2(token);
      T token_cast;
      ss2 >> token_cast;
      vec.push_back(token_cast);
  }
  return vec;
}

template<typename T>
std::ostream& operator<<(
std::ostream& os, const FixedWidthValue<T> value
)
{
  // Convert to string, using the same fmtflags as the ostream
  const auto os_flags = os.flags();
  std::stringstream ss;
  ss.flags(os_flags);
  ss << value.getValue();
  std::string raw_str( ss.str() );

  // Note: std::setw is not "sticky"
  const int width = value.getWidth();
  if ( width >= 0 ) {
    os << std::setw(width);
  }

  // Keep only part of string
  int len = raw_str.size();
  if ( width >= 0 and len > width ) {
    // Get 'adjustfield' flag, which determines text alignment
    std::ios_base::fmtflags adjust_field_flag = (os_flags & std::ios_base::adjustfield);

    // First character to keep
    int first = 0;
    if ( adjust_field_flag == std::ios_base::left ) {  // left-aligned
      first = 0;
    }
    else if ( adjust_field_flag == std::ios_base::right ) {
      first = len - width;
    }
    else {
      // TODO: how to consistently implement this?
      FANCY_ASSERT( adjust_field_flag == std::ios_base::internal, "unsupported alignment: internal");
      FANCY_ASSERT( false, "unrecognized adjustfield flag: " << adjust_field_flag );
    }

    os << raw_str.substr(first, width);
  }
  else {
    os << raw_str;
  }

  return os;
};

void swap(std::string& s1, std::string& s2);
} // end namespace StringTools

#endif // ifndef STRING_TOOLS_H
