#include "interface/interface.hpp"
#include "tools/Assert.hpp"
#include "tools/InputParser.hpp"
#include "analysis/factory.hpp"
int main(int argc, char **argv)
{
  FANCY_ASSERT(argc == 2, "Analysis code only accepts a single input that specifies the op input file.");
  std::string op_input_file_ = argv[1];
  std::string trajectory_file_, index_file_, topology_file_, gro_file_;
  InputParser input_parser;
  ParameterPack master_pack = input_parser.parseFile(op_input_file_);
  using KeyType = ParameterPack::KeyType;
  bool trajectory_found = master_pack.readString("trajectory", KeyType::Required, trajectory_file_);
  bool idx_found = master_pack.readString("index", KeyType::Optional, index_file_);
  bool top_found = master_pack.readString("top", KeyType::Optional, topology_file_);
  bool gro_found = master_pack.readString("gro", KeyType::Optional, gro_file_);
  FANCY_ASSERT(trajectory_found, "Failed to find trajectory file.");
  Box b1;
  if(idx_found) readNDX(index_file_, b1);
  if(top_found) readTOP(topology_file_, b1);
  if(gro_found) readGRO(gro_file_, b1);
  
  InputPack master_input_pack = InputPack(&master_pack, &b1);

  std::vector<InputPack> ag_packs = master_input_pack.buildDerivedInputPacks("AtomGroup");
  FANCY_ASSERT(ag_packs.size() > 0, "No atom groups specified.");
  for(std::size_t i = 0; i < ag_packs.size(); i++){
    std::string type, name;
    ag_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    ag_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto ag_ptr = AtomGroup_Factory(type, ag_packs[i]);
    master_input_pack.addAtomGroup(name, ag_ptr);
  } 

  std::vector<InputPack> pv_packs = master_input_pack.buildDerivedInputPacks("ProbeVolume");
  for(std::size_t i = 0; i < pv_packs.size(); i++){
    std::string type, name;
    pv_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    pv_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto pv_ptr = ProbeVolume_Factory(type, pv_packs[i]);
    master_input_pack.addProbeVolume(name, pv_ptr);
  }  
  
  std::vector<InputPack> calc_packs = master_input_pack.buildDerivedInputPacks("Calculation");
  for(std::size_t i = 0; i < calc_packs.size(); i++){
    std::string type, name;
    calc_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    calc_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto calc_ptr = Calculation_Factory(type, calc_packs[i]);
    master_input_pack.addCalculation(name, calc_ptr);
  }  

  std::cout << "Starting test..." << std::endl;
  XDRTrajectory traj(trajectory_file_);
  std::cout << "Trajectory loaded..." << std::endl;
  
  while(traj.nextFrame()){
    traj.getFrame(b1);
    std::cout << b1.time << "\n";
    auto ag_reg_ = master_input_pack.AtomGroupMap();
    auto pv_reg_ = master_input_pack.ProbeVolumeMap();
    auto calc_reg_ = master_input_pack.CalculationMap();
    //update all objects
    for(auto i = ag_reg_.begin(); i != ag_reg_.end(); i++){
      i->second->update();
    }
    for(auto i = pv_reg_.begin(); i != pv_reg_.end(); i++){
      i->second->update();
    }
    for(auto i = calc_reg_.begin(); i != calc_reg_.end(); i++){
      i->second->update();
    }
    //run all calculations
    for(auto i = calc_reg_.begin(); i != calc_reg_.end(); i++){
      i->second->calculate();
    }
    //perform all outputs
     for(auto i = calc_reg_.begin(); i != calc_reg_.end(); i++){
      i->second->printConsoleReport();
      i->second->output();
    }  


  }
  auto calc_reg_ = master_input_pack.CalculationMap();
  for(auto i = calc_reg_.begin(); i != calc_reg_.end(); i++){
    auto calculation = i->second;
    i->second->finalOutput();
  }

  return 0;
}