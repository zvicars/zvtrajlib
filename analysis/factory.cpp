#include "factory.hpp"
#include "../tools/Assert.hpp"

//PROBE VOLUMES
#include "probevolumes/PV_Simple.hpp"
#include "probevolumes/PV_Cylinder.hpp"
#include "probevolumes/PV_Boolean.hpp"
#include "probevolumes/PV_DynamicBox.hpp"
//ATOM GROUPS
#include "atomgroups/AtomGroup_name.hpp"
#include "atomgroups/AtomGroup_range.hpp"
#include "atomgroups/AtomGroup_resname.hpp"
#include "atomgroups/AtomGroup_ndx.hpp"
#include "atomgroups/AtomGroup_dynamic.hpp"
#include "atomgroups/AtomGroup_intersection.hpp"
#include "atomgroups/AtomGroup_union.hpp"
#include "atomgroups/AtomGroup_not.hpp"
#include "atomgroups/AtomGroup_system.hpp"
//CALCULATIONS
#include "calculations/Calc_Isosurface.hpp"
#include "calculations/Calc_Angle.hpp"
#include "calculations/Calc_Nv.hpp"
#include "calculations/Calc_Nv_wFit.hpp"
#include "calculations/Calc_Nv_wFit_IP.hpp"
#include "calculations/Calc_Dipole.hpp"
#include "calculations/Calc_2D_Density.hpp"
#include "calculations/Calc_1D_Density.hpp"
#include "calculations/Calc_1D_Density_IP.hpp"
#include "calculations/Calc_Write_Xtc.hpp"
#include "calculations/Calc_Lattice_Motion.hpp"
#include "calculations/Calc_Write_AvgPos.hpp"
#include "calculations/Calc_Relative_Pos.hpp"
#include "calculations/Calc_DensityField.hpp"
#include "calculations/Calc_DensityFieldExtra.hpp"
#include "calculations/Calc_DensityField_Electric.hpp"
#include "calculations/Calc_DensityField_Angle.hpp"
#include "calculations/Calc_IceID.hpp"

ProbeVolume* ProbeVolume_Factory(std::string key, InputPack& input){
  if(key == "rectilinear") return new PV_DiscreteRect(input);
  if(key == "dynbox") return new PV_DynBox(input);
  if(key == "cylinder") return new PV_Cylinder(input);
  if(key == "boolean") return new PV_Boolean(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}

AtomGroup* AtomGroup_Factory(std::string key, InputPack& input){
  if(key == "name") return new AtomGroup_name(input);
  if(key == "range") return new AtomGroup_range(input);
  if(key == "resname") return new AtomGroup_resname(input);
  if(key == "ndx") return new AtomGroup_ndx(input);
  if(key == "dynamic") return new AtomGroup_dynamic(input);
  if(key == "intersection") return new AtomGroup_intersection(input);
  if(key == "union") return new AtomGroup_union(input);
  if(key == "not") return new AtomGroup_not(input);
  if(key == "system") return new AtomGroup_system(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}

Calculation* Calculation_Factory(std::string key, InputPack& input){
  if(key == "isosurface") return new Calc_Isosurface(input);
  if(key == "nv") return new Calc_Nv(input);
  if(key == "nvwfit") return new Calc_Nv_wFit(input);
  if(key == "nvwfit_ip") return new Calc_Nv_wFit_IP(input);
  if(key == "dipole") return new Calc_Dipole(input);
  if(key == "angle") return new Calc_Angle(input);
  if(key == "3d_densityelectric") return new Calc_DensityFieldElectric(input);
  if(key == "3d_densityangle") return new Calc_DensityFieldAngle(input);
  if(key == "3d_densityextra") return new Calc_DensityFieldExtra(input);  
  if(key == "3d_density") return new Calc_DensityField(input);
  if(key == "2d_density") return new Calc_2D_Density(input);
  if(key == "1d_density") return new Calc_1D_Density(input);
  if(key == "1d_density_ip") return new Calc_1D_Density_IP(input);
  if(key == "circlefit") return new Calc_SWIPES_CircleFit(input);
  if(key == "writextc") return new Calc_Write_Xtc(input);
  if(key == "lattice_motion") return new Calc_Lattice_Motion(input);
  if(key == "write_avg_pos") return new Calc_Write_AvgPos(input);
  if(key == "rel_pos") return new Calc_Relative_Pos(input);
  if(key == "iceid") return new Calc_IceID(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}