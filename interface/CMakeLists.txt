add_subdirectory(libxdr)
add_subdirectory(ndx)
add_subdirectory(top)
add_library (TRAJLIB SHARED
datatypes.hpp
interface.hpp
interface.cpp
fileparsers/parseNDX.hpp
fileparsers/parseNDX.cpp
)
target_include_directories(TRAJLIB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TRAJLIB PRIVATE XDRFILE)
target_compile_options(TRAJLIB PRIVATE)
