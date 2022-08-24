add_executable(TestMatrixExe test_matrix.cpp)
target_link_libraries(TestMatrixExe PRIVATE TOOLS)
target_compile_options(TestMatrixExe PRIVATE -O0 -g)
add_test(NAME TestMatrix 
         COMMAND TestMatrixExe)