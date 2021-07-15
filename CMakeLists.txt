cmake_minimum_required(VERSION 3.0)
project(zvtrajlib VERSION 0.1)
set(CMAKE_CXX_STANDARD 11)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
add_subdirectory(libxdr)

add_executable(test test.cpp interface.cpp interface.h datatypes.h)
target_link_libraries(test PRIVATE XDRFILE)