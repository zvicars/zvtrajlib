#include "interface/interface.hpp"
#include "tools/Assert.hpp"
#include "tools/InputParser.hpp"
#include "analysis/InputPack.hpp"
#include "analysis/ProbeVolume.hpp"
#include "analysis/Calculation.hpp"

int main()
{
  const auto& calculate_registry = CalculationRegistry::Factory::factory().get_registry();
  if ( calculate_registry.size() < 1 ) {
    throw std::runtime_error("Error: Calculation registry is empty");
  }

  const auto& probe_volume_registry = ProbeVolumeRegistry::Factory::factory().get_registry();
  if ( probe_volume_registry.size() < 1 ) {
    throw std::runtime_error("Error: ProbeVolume registry is empty");
  }

  std::cout << "1" << std::endl;
  std::string op_input_file_ = "test.input";
  std::string trajectory_file_;
  InputParser input_parser;
  std::cout << "2" << std::endl;
  ParameterPack master_pack = input_parser.parseFile(op_input_file_);
  using KeyType = ParameterPack::KeyType;
  master_pack.readString("trajectory", KeyType::Required, trajectory_file_);
  InputPack master_input_pack = InputPack(&master_pack);
  std::cout << "3" << std::endl;
  std::vector<InputPack> pv_packs = master_input_pack.buildDerivedInputPacks("ProbeVolume");
  std::cout << "4" << std::endl;
  for(std::size_t i = 0; i < pv_packs.size(); i++){
    std::string type, name;
    pv_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    pv_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    std::cout << name << " " << type << std::endl;
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
  std::cout << "5" << std::endl;
  std::cout << "Starting test..." << std::endl;
  XDRTrajectory traj(trajectory_file_);
  Box b1;
  while(traj.nextFrame()){
    traj.getFrame(b1);
    auto calc_reg_ = master_input_pack.getCalculationRegistry();
    for(auto i = calc_reg_->begin(); i != calc_reg_->end(); i++)
      i->second->calculate(b1);
    }
  return 0;
}