add_subdirectory(golosio)
add_subdirectory(rchandra)
add_library (INSTANT_INT STATIC
)
target_sources(INSTANT_INT PRIVATE
MarchingCubesInterface.cpp
VoxelGrid.cpp
)
target_include_directories(INSTANT_INT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(INSTANT_INT PRIVATE GOLOSIO RCHANDRA)
target_compile_options(INSTANT_INT PRIVATE -O0)