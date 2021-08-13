add_subdirectory(golosio)
add_library (INSTANT_INT STATIC
MarchingCubesInterface.cpp
meshcurvature.cpp
VoxelGrid.cpp
)
target_include_directories(INSTANT_INT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../extern/eigen-3.3.9/)
target_link_libraries(INSTANT_INT PRIVATE GOLOSIO)
target_compile_options(INSTANT_INT PRIVATE -O0)

add_subdirectory(test)