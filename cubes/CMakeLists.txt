add_subdirectory(golosio)
add_library (INSTANT_INT STATIC
MarchingCubesInterface.cpp
meshcurvature.cpp
VoxelGrid.cpp
meshfiles/objfile.cpp
)


target_include_directories(INSTANT_INT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../extern/eigen-3.3.9/)
if(OpenMP_CXX_FOUND)
target_link_libraries(INSTANT_INT PRIVATE GOLOSIO OpenMP::OpenMP_CXX)
endif()
target_compile_options(INSTANT_INT PRIVATE -O0 -fopenmp)

add_subdirectory(test)