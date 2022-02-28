add_library (ANALYSIS STATIC "")
target_sources(ANALYSIS PRIVATE
bonded_functions.cpp
atomgroups/AtomGroup.cpp
atomgroups/AtomGroup_ndx.cpp
atomgroups/AtomGroup_name.cpp
atomgroups/AtomGroup_range.cpp
atomgroups/AtomGroup_resname.cpp
atomgroups/AtomGroup_dynamic.cpp
atomgroups/AtomGroup_intersection.cpp
atomgroups/AtomGroup_union.cpp
atomgroups/AtomGroup_not.cpp
atomgroups/AtomGroup_system.cpp
calculations/Calculation.cpp
calculations/Calc_Histogram.cpp
calculations/Calc_Nv.cpp
calculations/Calc_Nv_wFit.cpp
calculations/Calc_Nv_wFit_IP.cpp
calculations/Calc_Isosurface.cpp
calculations/Calc_Angle.cpp
calculations/Calc_DensityField.cpp
calculations/Calc_2D_Density.cpp
calculations/Calc_1D_Density.cpp
calculations/Calc_1D_Density_IP.cpp
calculations/Calc_SWIPES_CircleFit.cpp
calculations/Calc_Write_Xtc.cpp
calculations/Calc_Dipole.cpp
calculations/Calc_Lattice_Motion.cpp
calculations/Calc_Write_AvgPos.cpp
calculations/Calc_Relative_Pos.cpp
calculations/Calc_IceID.cpp
probevolumes/ProbeVolume.cpp
probevolumes/PV_Simple.cpp
probevolumes/PV_Cylinder.cpp
probevolumes/PV_Boolean.cpp
probevolumes/PV_DynamicBox.cpp
InputPack.cpp
factory.cpp
)
target_include_directories(ANALYSIS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../extern/eigen-3.3.9/)
target_link_libraries(ANALYSIS PRIVATE TRAJLIB
                               PRIVATE TOOLS
                               PRIVATE INSTANT_INT
                               PRIVATE XDRFILE)
if(OpenMP_CXX_FOUND)
target_link_libraries(ANALYSIS PRIVATE GOLOSIO OpenMP::OpenMP_CXX)
endif()
target_compile_options(ANALYSIS PRIVATE -O3 -fopenmp)