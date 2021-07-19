add_library (GOLOSIO STATIC "")
target_sources(GOLOSIO PRIVATE
TriangulateMain.c
Triangulate.c
TriangulateGlobal.c
FindVertexes.c
OutputData.c
PolygonizeCube.c
matrix.c
)
target_include_directories(GOLOSIO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(GOLOSIO PRIVATE -O3)