#include "Calc_HBonds.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/cellgrid.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../interface/interface.hpp"
inline double f_angle(Vec3<double> v21, Vec3<double> v23){
  double num = 0, norm1 = 0, norm2 = 0;
  for(int i = 0; i < 3; i++){
      num += v21[i]*v23[i];
      norm1 += v21[i]*v21[i];
      norm2 += v23[i]*v23[i];
  }
  double denom = sqrt(norm1)*sqrt(norm2);
  return acos(num/denom);
}

Calc_HBonds::Calc_HBonds(InputPack& input) : Calculation{input}{
  std::string atom_group_name, atom_group_name2;
  std::vector<std::string> h_names;
  input.params().readString("atom_group", KeyType::Required, atom_group_name);
  atom_group_ = input.findAtomGroup(atom_group_name);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atomgroup in Calc_Hbonds : " + name_);
  input.params().readNumber("r", KeyType::Required, r_);
  input.params().readNumber("angle", KeyType::Required, angle_);
  //length of vector decides the size of the water, 1 is the location of a hydrogen, 2 is the index of the + charge
  //tip4p example: water_format = [ 0 1 1 2 ]
  //spce example: water_format = [ 2 1 1 ]
  //computation will be done only on h-bond donor molecules
  std::vector<int> water_format; 
  input.params().readVector("water_format", KeyType::Required, water_format);
  for(int i = 0; i < water_format.size(); i++){
    if(water_format[i] == 2) rel_charge_index_ = i;
    if(water_format[i] == 1) rel_h_indices_.push_back(i);
  }
  input.params().readFlag("write_traj", KeyType::Required, write_traj_);
  if(write_traj_){
    std::ofstream ofile(name_+".xyz");
    ofile.close();
  }
  return;
}
void Calc_HBonds::calculate(){
  if(!doCalculate()) return;
  calculate_flag_ = 1;
  hbonds_.clear();
  auto& o_indices = atom_group_->getIndices();
  CellGrid c1(r_, box_size_);
  auto& atoms = box->atoms;
  #pragma omp parallel for
  for(auto index : o_indices){
    auto& atom = atoms[index];
    c1.addIndexToGrid(index, atom.x);
  }
  #pragma omp parallel
  {
    std::vector<HBond> hbond_private;
    #pragma omp for nowait
    for(auto index : o_indices){
      auto& atom = atoms[index];
      auto nearby_indices = c1.getNearbyIndices(index, atom.x);
      std::vector<Vec3<double> > donor_h_positions;
      std::vector<int> donor_h_indices;
      int donor_charge_index = index + rel_charge_index_;
      Vec3<double> x_od = atoms[donor_charge_index].x;
      for(auto rel_index : rel_h_indices_){
        donor_h_indices.push_back(index + rel_index);
        donor_h_positions.push_back(atoms[index + rel_index].x);
      }
      //looking for nearby h-atoms
      for(auto nearby_index : nearby_indices){
        auto x_oa = atoms[nearby_index + rel_charge_index_].x;
        int counter = -1; 
        for(auto x_hd : donor_h_positions){
          counter++;
          getNearestImage3D(x_oa, x_hd, box_size_);
          auto v_hdoa = x_oa-x_hd;
          auto v_hdod = x_od-x_hd;
          auto r = norm2(v_hdoa);
          //angle OD-HD-OA
          auto angle = f_angle(v_hdoa, v_hdod)*180/M_PI;
          if(r > r_) continue;
          if(angle < 180 - angle_) continue;
          HBond hbtemp;
          hbtemp.x = x_hd + (0.5*r*v_hdoa); hbtemp.v = v_hdoa; hbtemp.angle = angle;
          hbtemp.indices = {donor_h_indices[counter], nearby_index + rel_charge_index_}; hbtemp.r = r;
          hbond_private.push_back(hbtemp);
        }
      }
    }
    #pragma omp critical
    hbonds_.insert(hbonds_.end(), hbond_private.begin(), hbond_private.end());
  }
  return;
}
void Calc_HBonds::finalOutput(){
  return;
}
void Calc_HBonds::update(){
  Calculation::update();
  auto bv = box->boxvec;
  for(int i = 0; i < 3; i++){
    box_size_[i] = bv[i][i];
  }
  return;
}
void Calc_HBonds::output(){
  //need to use xyz format for this
  if(!doOutput()) return;
  if(!write_traj_) return; 
  std::vector<std::string> names;
  names.reserve(hbonds_.size()*2);
  std::vector<Vec3<double> > positions;
  positions.reserve(hbonds_.size()*2);
  for(auto hb : hbonds_){
    positions.push_back(hb.x);
    positions.push_back(hb.x + hb.v);
    names.push_back("T");
    names.push_back("H");
  }
  writeXYZ_ov_append(name_+".xyz", box_size_, names, positions);
  return;
}