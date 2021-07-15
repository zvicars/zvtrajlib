cmake_minimum_required(VERSION 3.0)
project(zvtrajlib VERSION 0.1)
set(CMAKE_CXX_STANDARD 11)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
add_subdirectory(libxdr)
add_library (ZVTRAJLIB 
bonded_functions.h
datatypes.h
interface.h
interface.cpp
)
target_include_directories(ZVTRAJLIB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ZVTRAJLIB PRIVATE XDRFILE)
add_executable(test test.cpp)
target_link_libraries(test PRIVATE ZVTRAJLIB)
