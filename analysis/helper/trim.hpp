#include <iostream>
#include <string>
#include <algorithm>
 
const std::string WHITESPACE = " \n\r\t\f\v";
 
static inline std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
static inline std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
 
static inline std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}