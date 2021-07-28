#pragma once
#include "../tools/Assert.hpp"
#include "probevolumes/PV_Simple.hpp"
#include "atomgroups/AtomGroup_name.hpp"
#include "atomgroups/AtomGroup_ndx.hpp"
#include "calculations/Calc_Isosurface.hpp"
#include "calculations/Calc_Angle.hpp"
#include "calculations/Calc_Nv.hpp"

ProbeVolume* ProbeVolume_Factory(std::string key, InputPack& input);
AtomGroup* AtomGroup_Factory(std::string key, InputPack& input);
Calculation* Calculation_Factory(std::string key, InputPack& input);