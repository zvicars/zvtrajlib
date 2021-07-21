#include "interface/interface.hpp"
#include "tools/Assert.hpp"
#include "tools/InputParser.hpp"
#include "analysis/InputPack.hpp"
#include "analysis/ProbeVolume.hpp"
#include "analysis/Calculation.hpp"
#include "analysis/AtomGroup.hpp"
void checkRegistries(){
  const auto& probe_volume_registry = ProbeVolumeRegistry::Factory::factory().get_registry();
  if ( probe_volume_registry.size() < 1 ) {
    throw std::runtime_error("Error: ProbeVolume registry is empty");
  }
  const auto& calculate_registry = CalculationRegistry::Factory::factory().get_registry();
  if ( calculate_registry.size() < 1 ) {
    throw std::runtime_error("Error: Calculation registry is empty");
  }
  const auto& atomgroup_registry = AtomGroupRegistry::Factory::factory().get_registry();
  if ( calculate_registry.size() < 1 ) {
    throw std::runtime_error("Error: AtomGroup registry is empty");
  }
  return;
}

int main(int argc, char **argv)
{
  checkRegistries();
  FANCY_ASSERT(argc == 2, "Analysis code only accepts a single input that specifies the op input file.");
  std::string op_input_file_ = argv[1];
  std::string trajectory_file_, index_file_, topology_file_;
  InputParser input_parser;
  ParameterPack master_pack = input_parser.parseFile(op_input_file_);
  using KeyType = ParameterPack::KeyType;
  bool trajectory_found = master_pack.readString("trajectory", KeyType::Required, trajectory_file_);
  bool idx_found = master_pack.readString("index", KeyType::Optional, index_file_);
  bool top_found = master_pack.readString("top", KeyType::Optional, topology_file_);
  FANCY_ASSERT(trajectory_found, "Failed to find trajectory file.");
  Box b1;
  if(idx_found) readNDX(index_file_, b1);
  if(top_found) readTOP(topology_file_, b1);
  InputPack master_input_pack = InputPack(&master_pack, &b1);

  std::vector<InputPack> ag_packs = master_input_pack.buildDerivedInputPacks("AtomGroup");
  FANCY_ASSERT(ag_packs.size() > 0, "No atom groups specified.");
  for(std::size_t i = 0; i < ag_packs.size(); i++){
    std::string type, name;
    ag_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    ag_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto ag_ptr = AtomGroupRegistry::Factory::factory().create(type, ag_packs[i]);
    master_input_pack.addAtomGroup(name, ag_ptr);
  } 

  std::vector<InputPack> pv_packs = master_input_pack.buildDerivedInputPacks("ProbeVolume");
  for(std::size_t i = 0; i < pv_packs.size(); i++){
    std::string type, name;
    pv_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    pv_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto pv_ptr = ProbeVolumeRegistry::Factory::factory().create(type, pv_packs[i]);
    master_input_pack.addProbeVolume(name, pv_ptr);
  }  
  
  std::vector<InputPack> calc_packs = master_input_pack.buildDerivedInputPacks("Calculation");
  for(std::size_t i = 0; i < calc_packs.size(); i++){
    std::string type, name;
    calc_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    calc_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto calc_ptr = CalculationRegistry::Factory::factory().create(type, calc_packs[i]);
    master_input_pack.addCalculation(name, calc_ptr);
  }  

  std::cout << "Starting test..." << std::endl;
  XDRTrajectory traj(trajectory_file_);
  std::cout << "Trajectory loaded..." << std::endl;
  
  while(traj.nextFrame()){
    traj.getFrame(b1);
    auto calc_reg_ = master_input_pack.CalculationMap();
    for(auto i = calc_reg_.begin(); i != calc_reg_.end(); i++){
      auto calculation = i->second;
      i->second->calculate();
      //std::cout << i->second->printConsoleReport();
    }
  }
  return 0;
}