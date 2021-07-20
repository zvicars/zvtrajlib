add_library(STUFF STATIC
marching_cubes/Triangulate.c
marching_cubes/TriangulateGlobal.c
marching_cubes/FindVertexes.c
marching_cubes/OutputData.c
marching_cubes/PolygonizeCube.c
marching_cubes/matrix.c
marching_cubes/matrix.h
marching_cubes/Triangulate.h
marching_cubes/TriangulateGlobal.h
marching_cubes/TriangulateTab.h
xdr_interface.h
xdr_wrapper.h
marching_cubes.h
parse_input.h
basic_structs.h
../../libxdr/trr_seek.h
../../libxdr/trr_seek.c
../../libxdr/xdrfile_trr.h
../../libxdr/xdrfile_trr.c
../../libxdr/xdrfile_xtc.h
../../libxdr/xdrfile_xtc.c
../../libxdr/xdrfile.h
../../libxdr/xdrfile.c
../../libxdr/xtc_seek.h
../../libxdr/xtc_seek.c

)
target_include_directories(STUFF PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/libxdr)

add_executable(old main.cpp)
target_link_libraries(old PRIVATE STUFF)
target_compile_options(old PRIVATE -O0)