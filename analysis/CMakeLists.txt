add_library (ANALYSIS
bonded_functions.hpp
bonded_functions.cpp
calculations.hpp
calculations.cpp
ProbeVolume.hpp
ProbeVolume.cpp
)
target_include_directories(ANALYSIS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ANALYSIS PRIVATE INTERFACE IO)