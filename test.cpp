#include "interface/interface.hpp"
#include "analysis/ProbeVolume.hpp"
#include "analysis/Calculation.hpp"
int main()
{
  std::string op_input_file_ = "test.input";
  std::string trajectory_file_;
  InputParser input_parser;

  ParameterPack master_pack = input_parser.parseFile(op_input_file_);
  using KeyType = ParameterPack::KeyType;
  master_pack.readString("trajectory", KeyType::Required, trajectory_file_);
  InputPack master_input_pack = InputPack(&master_pack);

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
  Box b1;
  while(traj.nextFrame()){
    traj.getFrame(b1);
    auto calc_reg_ = master_input_pack.getCalculationRegistry();
    for(auto i = calc_reg_->begin(); i != calc_reg_->end(); i++)
      i->second->calculate(b1);
    }
  return 0;
}