add_library (TOOLS STATIC
  Stacktrace.cpp
  InputParser.cpp
  StringTools.cpp
  meshtools_files/meshops.cpp
  meshtools_files/obj.cpp
)
target_include_directories(TOOLS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(TOOLS PRIVATE -fPIC)

add_subdirectory(test)