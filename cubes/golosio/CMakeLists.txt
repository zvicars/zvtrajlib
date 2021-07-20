add_library (GOLOSIO STATIC 
Triangulate.cpp
TriangulateGlobal.cpp
FindVertexes.cpp
PolygonizeCube.cpp
matrix.cpp
TriangulateGlobal.hpp
)
target_include_directories(GOLOSIO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(GOLOSIO PRIVATE -O0)