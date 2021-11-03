#include "factory.hpp"
#include "../tools/Assert.hpp"
#include "probevolumes/PV_Simple.hpp"
#include "probevolumes/PV_Cylinder.hpp"
#include "probevolumes/PV_Boolean.hpp"
#include "atomgroups/AtomGroup_name.hpp"
#include "atomgroups/AtomGroup_resname.hpp"
#include "atomgroups/AtomGroup_ndx.hpp"
#include "atomgroups/AtomGroup_dynamic.hpp"
#include "atomgroups/AtomGroup_intersection.hpp"
#include "atomgroups/AtomGroup_union.hpp"
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
ProbeVolume* ProbeVolume_Factory(std::string key, InputPack& input){
  if(key == "rectilinear") return new PV_DiscreteRect(input);
  if(key == "cylinder") return new PV_Cylinder(input);
  if(key == "boolean") return new PV_Boolean(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}

AtomGroup* AtomGroup_Factory(std::string key, InputPack& input){
  if(key == "name") return new AtomGroup_name(input);
  if(key == "resname") return new AtomGroup_resname(input);
  if(key == "ndx") return new AtomGroup_ndx(input);
  if(key == "dynamic") return new AtomGroup_dynamic(input);
  if(key == "intersection") return new AtomGroup_intersection(input);
  if(key == "union") return new AtomGroup_union(input);
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
  if(key == "2d_density") return new Calc_2D_Density(input);
  if(key == "1d_density") return new Calc_1D_Density(input);
  if(key == "1d_density_ip") return new Calc_1D_Density_IP(input);
  if(key == "circlefit") return new Calc_SWIPES_CircleFit(input);
  if(key == "writextc") return new Calc_Write_Xtc(input);
  if(key == "lattice_motion") return new Calc_Lattice_Motion(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}