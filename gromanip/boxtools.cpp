#include "boxtools.hpp"
#include <set>
#include <algorithm>
#include <unordered_set>
#include "../tools/StringTools.hpp"
#include "Eigen/Eigen"
bool atom_res_sorter(const Atom& lhs, const Atom& rhs){
  if(lhs.resnr != rhs.resnr) return lhs.resnr < rhs.resnr;
  return lhs.index < rhs.index;
}

bool boxtools::renumberBox(Box& box){
  boxtools::checkBoxValid(box);
  std::sort(box.atoms.begin(), box.atoms.end(), &atom_res_sorter);
  int prevResNum = box.atoms[0].resnr;
  int resCounter = 1;
  int atomCounter = 1;
  for(auto & atom : box.atoms){
    if(atom.resnr != prevResNum){
      resCounter++;
      prevResNum = atom.resnr;
    }
    atom.resnr = resCounter;
    atom.index = atomCounter;
    atomCounter++;
  }
  return 0;
}

bool offsetRes(Box& box, const Atom& last){
  for(auto & atom : box.atoms){
    atom.resnr += last.resnr;
    atom.index += last.index;
  }
  return 0;
}
bool boxtools::deleteRestype(Box& box, std::string resname){
  boxtools::checkBoxValid(box);
  auto it = box.atoms.begin();
  while (it != box.atoms.end()){
    if (it->resname == resname){
        it = box.atoms.erase(it);
    }
    else {
        ++it;
    }    
  }
  return 0;
}
Box boxtools::mergeBox(const Box& b1, const Box& b2){
  boxtools::checkBoxValid(b1);
  boxtools::checkBoxValid(b2);
  Box retBox = b1;
  Box b2_temp = b2;
  //should ensure that resNum of b2 is always > resNum of b1, will ensure proper numbering at end
  offsetRes(b2_temp, b1.atoms.back());
  for(int i = 0; i < 3; i++){
    retBox.boxvec[i][i] = std::max(b1.boxvec[i][i], b2.boxvec[i][i]);
  }
  retBox.atoms.insert(retBox.atoms.begin(), b2_temp.atoms.begin(), b2_temp.atoms.end());
  renumberBox(retBox);
  return retBox;
}

void boxtools::removeResNumbers(Box& box, std::vector<int> res_list){
  boxtools::checkBoxValid(box);
  //make sure all resnumbers are unique
  std::unordered_set<int> res_list_set;
  for (const int &i: res_list) {
    res_list_set.emplace(i);
  } 
  std::vector<Atom> new_atoms;
  new_atoms.reserve(box.atoms.size());
  //remove all resnumbers from res_list_set
  for(auto& atom : box.atoms){
    bool remove_atom = 0;
    for(int resnum : res_list_set){
      if(resnum == atom.resnr){
        remove_atom = 1;
        break;
      }
    }
    if(!remove_atom) new_atoms.push_back(atom);
  }
  box.atoms = new_atoms;
  renumberBox(box);
  return;
}

std::vector<int> boxtools::getResnrWithinVolume(const Box& box, const Volume& volume){
  std::vector<int> retVec;
  for(auto& atom : box.atoms){
    if(!volume.isInside(atom.x)) retVec.push_back(atom.resnr);
  }
  return retVec;
}
std::vector<int> boxtools::getResnrNotWithinVolumebyAtomName(const Box& box, const Volume& volume, std::string at_name){
  std::vector<int> retVec;
  for(auto& atom : box.atoms){
    if(!volume.isInside(atom.x) && atom.name == at_name) retVec.push_back(atom.resnr);
  }
  return retVec;
}
std::vector<int> boxtools::getResnrWithinVolumebyAtomName(const Box& box, const Volume& volume, std::string at_name){
  std::vector<int> retVec;
  for(auto& atom : box.atoms){
    if(volume.isInside(atom.x) && atom.name == at_name) retVec.push_back(atom.resnr);
  }
  return retVec;
}
std::vector<int> boxtools::getResnrAllButResname(const Box& box, const Volume& volume, std::string res_name){
  std::vector<int> retVec;
  for(auto& atom : box.atoms){
    if(!volume.isInside(atom.x) && atom.resname == res_name) retVec.push_back(atom.resnr);
  }
  return retVec;
}

void boxtools::rotateEulerAngles(Box& box, Vec3<double> angles){
  Eigen::Vector3d radangles;
  radangles << angles[2]*M_PI/180.0, angles[1] * M_PI / 180.0, angles[0] * M_PI / 180.0;
  Eigen::AngleAxisd rollAngle(radangles[0], Eigen::Vector3d::UnitZ());
  Eigen::AngleAxisd yawAngle(radangles[1], Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd pitchAngle(radangles[2], Eigen::Vector3d::UnitX());
  Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;
  Eigen::Matrix3d rotationMatrix = q.matrix();
  for(auto& atom : box.atoms){
    Eigen::Vector3d xi, xi2;
    xi << atom.x[0], atom.x[1], atom.x[2];
    xi2 = rotationMatrix*xi;
    atom.x[0] = xi2[0];
    atom.x[1] = xi2[1];
    atom.x[2] = xi2[2];
  }
  return;
}

void boxtools::translateAtoms(Box& box, Vec3<double> offset){
  for(auto& atom : box.atoms){
    atom.x[0] += offset[0];
    atom.x[1] += offset[1];
    atom.x[2] += offset[2];
  }
  return;
}

//move corner of bounding box to 0,0,0 and shrink box dimensions to fit atoms
void boxtools::shrinkWrap(Box& box){
  Vec3<double> min, max;
  min.fill(std::numeric_limits<double>::max());
  max.fill(std::numeric_limits<double>::min());
  for(auto& atom : box.atoms){
    for(int i = 0; i < 3; i++){
      if(atom.x[i] < min[i]) min[i] = atom.x[i];
      if(atom.x[i] > max[i]) max[i] = atom.x[i];
    }
  }
  Vec3<double> newMin = {0};
  Vec3<double> newMax;
  for(int i = 0; i < 3; i++){
    newMax[i] = max[i] - min[i];
  }

  for(auto& atom : box.atoms){
    for(int i = 0; i < 3; i++){
      atom.x[i] -= min[i];
    }
  }
  for(int i = 0; i < 3; i++){
    box.boxvec[i][i] = newMax[i];
  }
  renumberBox(box);
  return;
}

void boxtools::relabelRes(Box& box, std::string original_name, std::string new_name){
  Vec3<double> min = {std::numeric_limits<double>::max()};
  Vec3<double> max = {std::numeric_limits<double>::min()};
  for(auto& atom : box.atoms){
    if(atom.resname == original_name){
      atom.resname = new_name;
    }
  }
  return;
}

void boxtools::relabelAtom(Box& box, std::string original_name, std::string new_name){
  Vec3<double> min = {std::numeric_limits<double>::max()};
  Vec3<double> max = {std::numeric_limits<double>::min()};
  for(auto& atom : box.atoms){
    if(atom.name == original_name){
      atom.name = new_name;
    }
  }
  return;
}

void boxtools::setBoxSize(Box& box, Vec3<double> newsize){
  for(int i = 0; i < 3; i++){
    box.boxvec[i][i] = newsize[i];
  }
  return;
}

