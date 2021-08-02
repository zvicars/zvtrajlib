add_executable(TestVoxelGridExe voxelgridtest.cpp)
target_link_libraries(TestVoxelGridExe PRIVATE INSTANT_INT)
target_compile_options(TestVoxelGridExe PRIVATE -O0 -g)
add_test(NAME TestVoxelGrid 
         COMMAND TestVoxelGridExe)