add_executable(TestVoxelGridExe voxelgridtest.cpp)
target_link_libraries(TestVoxelGridExe PRIVATE INSTANT_INT)
target_compile_options(TestVoxelGridExe PRIVATE -O0 -g)
add_test(NAME TestVoxelGrid 
         COMMAND TestVoxelGridExe)

add_executable(TestCurvatureExe testcurvature.cpp)
target_link_libraries(TestCurvatureExe PRIVATE INSTANT_INT)
target_compile_options(TestCurvatureExe PRIVATE -O0 -g)
add_test(NAME TestCurvature
        COMMAND TestCurvatureExe)

file(COPY testfiles DESTINATION ${CMAKE_CURRENT_BINARY_DIR})