add_library (TOOLS 
InputParser.hpp
InputParser.cpp
StringTools.cpp
StringTools.hpp
Stacktrace.hpp
Assert.hpp
io_headers.hpp
string_ops.hpp
GenericFactory.hpp
)
target_include_directories(TOOLS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})