add_library (ANALYSIS OBJECT "")
target_sources(ANALYSIS PRIVATE
bonded_functions.cpp
Calculation.cpp
Calc_Nv.cpp
ProbeVolume.cpp
PV_Simple.cpp
InputPack.hpp
Calculation.hpp
Calc_Nv.hpp
ProbeVolume.hpp
PV_Simple.hpp
)
target_include_directories(ANALYSIS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ANALYSIS PRIVATE TRAJLIB
                               PRIVATE TOOLS
                               PRIVATE INSTANT_INT)
target_compile_options(ANALYSIS PRIVATE -fPIC)