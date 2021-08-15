add_library (GOLOSIO STATIC
TriangulateGlobal.hpp
TriangulateGlobal.cpp
FindVertexes.cpp
PolygonizeCube.cpp
matrix.cpp
)
target_include_directories(GOLOSIO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(GOLOSIO PRIVATE -O3)