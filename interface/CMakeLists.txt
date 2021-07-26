add_subdirectory(libxdr)
add_library (TRAJLIB SHARED
interface.cpp
fileparsers/parseNDX.cpp
fileparsers/parseGRO.cpp
)
target_include_directories(TRAJLIB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TRAJLIB PRIVATE XDRFILE)
target_compile_options(TRAJLIB PRIVATE -O0)
