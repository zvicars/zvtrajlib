add_library (IO 
InputParser.hpp
InputParser.cpp
StringTools.cpp
StringTools.hpp
Stacktrace.hpp
Assert.hpp
io_headers.hpp
string_ops.hpp
)
target_include_directories(IO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})