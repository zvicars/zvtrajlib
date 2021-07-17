add_library (IO
InputParser.hpp
InputParser.cpp
io_headers.h
string_ops.h
)
target_include_directories(IO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})