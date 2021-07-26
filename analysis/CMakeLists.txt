add_library (ANALYSIS OBJECT "")
target_sources(ANALYSIS PRIVATE
bonded_functions.cpp
AtomGroup.cpp
atomgroups/AtomGroup_ndx.cpp
atomgroups/AtomGroup_name.cpp
Calculation.cpp
calculations/Calc_Nv.cpp
calculations/Calc_Isosurface.cpp
calculations/Calc_Angle.cpp
ProbeVolume.cpp
probevolumes/PV_Simple.cpp
)
target_include_directories(ANALYSIS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ANALYSIS PRIVATE TRAJLIB
                               PRIVATE TOOLS
                               PRIVATE INSTANT_INT)
target_compile_options(ANALYSIS PRIVATE -O0)