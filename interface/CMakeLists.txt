add_library (ZVTRAJLIB 
datatypes.hpp
interface.hpp
interface.cpp
)
target_include_directories(ZVTRAJLIB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ZVTRAJLIB PRIVATE XDRFILE)

