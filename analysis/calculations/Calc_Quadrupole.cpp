#include "Calc_Quadrupole.hpp"
#include "Eigen/Eigen"
Calc_Quadrupole::Calc_Quadrupole(InputPack& input):Calculation{input}
{
  std::string agname;
  input.params().readString("atomgroup", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  Vec<double> normal;
  input.params().readVector("normal", KeyType::Required, normal); //surface normal for cos theta calc
  FANCY_ASSERT(normal.size() == 3, "Improperly sized normal vector in quadrupole calc");
  for(int i = 0; i < 3; i++){
    normal_[i] = normal[i];
  }
  FANCY_ASSERT(atom_group_ != 0, "Failed to find atom group.");
  FANCY_ASSERT(atom_group_->getType() == "resname", "Can only use resname type at the moment.");
  auto molinfo = input.params().findParameterPacks("molecule_info", KeyType::Required);
  FANCY_ASSERT(molinfo.size() == 1, "Expecting only a single molinfo object in the calculation");
  molinfo[0]->readNumber("natoms", KeyType::Required, natoms_);
  molinfo[0]->readVector("charges", KeyType::Required, charges_);
  molinfo[0]->readVector("masses", KeyType::Required,  masses_);
  FANCY_ASSERT(charges_.size() == natoms_, "Improperly sized charge vector");
  FANCY_ASSERT(masses_.size() == natoms_, "Improperly sized mass vector");
  return;
}

//quadrupole calculation will be done on a per-molecule basis
//step1: calculate com and proper periodic image for each atom in a molecule
//step2: construct a quadrupole tensor for each molecule and extract the largest eigenvector

Eigen::Vector3d get_largest_eigen_vector(Eigen::Matrix3d A){
  Eigen::EigenSolver<Eigen::Matrix3d> es(A);
  auto eigenvalue = es.eigenvectors().col(0);
  return eigenvalue;
}

void Calc_Quadrupole::calculate(){
  auto indices = atom_group_->getIndices();
  int totAtoms = indices.size();
  int mol_it = 0, at_it = 0, mol_at_it;
  Eigen::Matrix3d qtensor;
  Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
  std::vector<Eigen::Vector3d> eigen_vecs(totAtoms/natoms_);
  while(mol_it < totAtoms/natoms_)
  {
    auto pos = box->atoms[indices[at_it]].x;
    Eigen::Vector3d eig_pos;
    eig_pos << pos[0], pos[1], pos[2];
    double mag2 = eig_pos.squaredNorm();
    qtensor = qtensor + charges_[mol_at_it]*(3*eig_pos*eig_pos.transpose() - mag2*identity);
    //molecule is fully specified, extract relevant info
    if(at_it%natoms_ == natoms_-1){
      auto vec = get_largest_eigen_vector(qtensor);
      eigen_vecs[mol_it] = vec;
      mol_it++;
      mol_at_it = 0;
    }
    at_it++;
    mol_at_it++;
  }
  //should have calculated an axis for each molecule, now calculate the cos(theta) between that and the surface normal
  

}

void Calc_Quadrupole::update(){
  if(hasUpdated()) return;
  Calculation::update();
  return; 
}
std::string Calc_Quadrupole::printConsoleReport(){
  return;
}

void Calc_Quadrupole::finalOutput(){
  if(output_freq_ <= 0) return;
  return;
}