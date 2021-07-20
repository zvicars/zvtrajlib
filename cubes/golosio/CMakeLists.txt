add_library (GOLOSIO OBJECT "")
target_sources(GOLOSIO PRIVATE
Triangulate.c
TriangulateGlobal.c
FindVertexes.c
OutputData.c
PolygonizeCube.c
matrix.c
)
target_include_directories(GOLOSIO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(GOLOSIO PRIVATE -O4)