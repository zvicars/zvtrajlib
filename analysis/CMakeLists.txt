add_library (ANALYSIS
bonded_functions.hpp
bonded_functions.cpp
calculations.hpp
calculations.cpp
parser.hpp
parser.cpp
probevolumes.hpp
probevolumes.cpp
)
target_include_directories(IO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(IO PRIVATE INTERFACE IO)