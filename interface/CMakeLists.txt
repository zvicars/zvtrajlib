add_library (ZVTRAJLIB 
bonded_functions.h
datatypes.h
interface.h
interface.cpp
utility.h
)
target_include_directories(ZVTRAJLIB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ZVTRAJLIB PRIVATE XDRFILE)

