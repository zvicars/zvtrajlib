add_library (RCHANDRA STATIC
CIsoSurface.hpp
Vectors.cpp
CIsoSurface.cpp
)
target_include_directories(RCHANDRA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(RCHANDRA PRIVATE -O0)