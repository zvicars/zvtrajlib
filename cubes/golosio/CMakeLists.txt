add_library (GOLOSIO STATIC 
Triangulate.c
TriangulateGlobal.c
FindVertexes.c
OutputData.c
PolygonizeCube.c
matrix.c
matrix.h
Triangulate.h
TriangulateGlobal.h
TriangulateTab.h
)
target_include_directories(GOLOSIO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(GOLOSIO PRIVATE -O0)