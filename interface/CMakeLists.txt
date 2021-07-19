add_library (TRAJLIB SHARED
datatypes.hpp
interface.hpp
interface.cpp
)
target_include_directories(TRAJLIB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TRAJLIB PRIVATE XDRFILE)
target_compile_options(TRAJLIB PRIVATE -fPIC)
