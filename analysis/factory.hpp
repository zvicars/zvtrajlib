#pragma once
#include "atomgroups/AtomGroup.hpp"
#include "probevolumes/ProbeVolume.hpp"
#include "calculations/Calculation.hpp"
#include "calculations/Calc_SWIPES_CircleFit.hpp"

ProbeVolume* ProbeVolume_Factory(std::string key, InputPack& input);
AtomGroup* AtomGroup_Factory(std::string key, InputPack& input);
Calculation* Calculation_Factory(std::string key, InputPack& input);