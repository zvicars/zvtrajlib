add_subdirectory(golosio)
add_library (INSTANT_INT STATIC
MarchingCubesInterface.hpp
MarchingCubesInterface.cpp
VoxelGrid.cpp
)
target_include_directories(INSTANT_INT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(INSTANT_INT PRIVATE GOLOSIO)
target_compile_options(INSTANT_INT PRIVATE -O0)

add_subdirectory(test)