cmake_minimum_required(VERSION 3.0)
enable_testing()
project(zvtrajlib VERSION 0.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)
add_compile_options(-O0)
find_package(OpenMP)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
add_subdirectory(interface)
add_subdirectory(tools)
add_subdirectory(cubes)
add_subdirectory(analysis)
add_subdirectory(gromanip)
add_executable(zvtrajlib zvtrajlib.cpp)
add_executable(zvmangro zvmangro.cpp)

target_link_libraries(zvtrajlib PRIVATE TRAJLIB 
                           PRIVATE TOOLS 
                           PRIVATE ANALYSIS)
target_compile_options(zvtrajlib PRIVATE -O0 -g)

target_link_libraries(zvmangro PRIVATE GROMANIP)
target_compile_options(zvmangro PRIVATE -O0 -g)