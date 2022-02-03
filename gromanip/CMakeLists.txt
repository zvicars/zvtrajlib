add_library (GROMANIP STATIC
  boxtools.cpp
  boxtools_files.cpp
  make_params.cpp  
  actions.cpp
  actions_crystals.cpp
  GroManipData.cpp
)

target_include_directories(GROMANIP PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../extern/eigen-3.3.9/)
target_link_libraries(GROMANIP PRIVATE TOOLS TRAJLIB)
target_compile_options(GROMANIP PRIVATE -O3)