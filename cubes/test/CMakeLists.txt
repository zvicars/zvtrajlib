add_executable(TestVoxelGridExe voxelgridtest.cpp)
target_link_libraries(test PRIVATE INSTANT_INT)
target_compile_options(test PRIVATE -O0 -g)
add_test(NAME TestVoxelGrid 
         COMMAND TestVoxelGridExe)