add_library (ANALYSIS
bonded_functions.hpp
bonded_functions.cpp
Calculation.hpp
Calculation.cpp
ProbeVolume.hpp
ProbeVolume.cpp
PV_Simple.cpp
PV_Simple.hpp
InputPack.hpp
)
target_include_directories(ANALYSIS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ANALYSIS PRIVATE INTERFACE TOOLS)