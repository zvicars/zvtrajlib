#include "boxtools.hpp"
#include <set>
double getPeriodicDistance1D(double x1, double x2, double lx){
	double dx = fabs(x2-x1);
	//nested ternary operators, in my code!? it's more likely than you think
	double eval = (dx > 0.5*lx) ? (dx-lx) : dx;
	return eval;
}
double getPeriodicDistance2(const Atom& a1, const Atom& a2, const Vec3<double>& box_dims, const Vec3<bool>& periodic_directions){
	//structure of box_dims is v1x, v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y) where v1 v2 and v3 are the box vectors
	//oh boy I sure do love ternary operators
	double dx = periodic_directions[0] ? getPeriodicDistance1D(a1.x[0], a2.x[0], box_dims[0]) : (a2.x[0] - a1.x[0]);
	double dy = periodic_directions[1] ? getPeriodicDistance1D(a1.x[1], a2.x[1], box_dims[1]) : (a2.x[1] - a1.x[1]);
	double dz = periodic_directions[2] ? getPeriodicDistance1D(a1.x[2], a2.x[2], box_dims[2]) : (a2.x[2] - a1.x[2]);
	return dx*dx + dy*dy + dz*dz;
}
int isBondInTable(const Atom& i, const Atom& j, const boxtools::ParameterTable<BondType>& btable)
{
  int counter = 0;
  for(auto& bond : btable.entries){
    if((bond.i == i.name && bond.j == j.name) || (bond.i == j.name && bond.j == i.name)) return counter;
    counter++;
  }
  return -1;
}

int isConstraintInTable(const Atom& i, const Atom& j, const boxtools::ParameterTable<ConstraintType>& btable)
{
  int counter = 0;
  for(auto& bond : btable.entries){
    if((bond.i == i.name && bond.j == j.name) || (bond.i == j.name && bond.j == i.name)) return counter;
    counter++;
  }
  return -1;
}

bool isExclusionInTable(const Atom& i, const Atom& j, const boxtools::ParameterTable<ExclusionType>& t)
{
  int counter = 0;
  for(auto& exclusion : t.entries){
    for(auto& name : exclusion.exclusion_names){
      if(exclusion.i == i.name && name == j.name) return 1;
      if(exclusion.i == j.name && name == i.name) return 1;
    }
    counter++;
  }
  return 0;
}
int isAtomInTable(const Atom& i, const boxtools::ParameterTable<AtomType>& atable)
{
  int counter = 0;
  for(auto& atom : atable.entries){
    if(atom.name == i.name) return counter;
    counter++;
  }
  return -1;
}
int isAngleInTable(const Atom& i, const Atom& j, const Atom& k, const boxtools::ParameterTable<AngleType>& atable, int& idx)
{
  int counter = 0;
  for(auto angle : atable.entries){
    // 0 no swapping needed ijk 
    // 1 swap indices 2,3 ikj 
    // 2 swap indices 1,2 jik 
    if(i.name == angle.i && j.name == angle.j && k.name == angle.k){idx = counter; return 0;}
    if(i.name == angle.i && j.name == angle.k && k.name == angle.j){idx = counter; return 1;}
    if(i.name == angle.j && j.name == angle.i && k.name == angle.k){idx = counter; return 2;}
    if(i.name == angle.j && j.name == angle.k && k.name == angle.i){idx = counter; return 1;}
    if(i.name == angle.k && j.name == angle.i && k.name == angle.j){idx = counter; return 2;}
    if(i.name == angle.k && j.name == angle.j && k.name == angle.i){idx = counter; return 0;}
    counter++;
  }
  idx = -1;
  return -1;
}

//require an exact match due to an unconstrained search (remains to be seen how slow this is)
int isVsiteInTable(std::string a, std::string b, std::string c, std::string d, const boxtools::ParameterTable<Vsite3Type>& t){
  int pos = 0;
  std::vector<std::string> names;
  names.push_back(b); names.push_back(c); names.push_back(d);  
  for(auto vsite : t.entries){
    if(vsite.i != a) continue;
    std::vector<std::string> ref_names;
    ref_names.push_back(vsite.j); ref_names.push_back(vsite.k); ref_names.push_back(vsite.l);
    if(names == ref_names) return pos;
    pos++;
  }  
  return -1;
}

std::vector<AtomInst> boxtools::makeAtoms(const Box& box, const boxtools::ParameterTable<BTAtomType>& t){
  std::vector<AtomInst> output;
  if( t.entries.size() == 0) return output;
  for(const auto& atom : box.atoms){
    for(auto& atomtype : t.entries){
      if(atom.name != atomtype.name) continue;
      AtomInst new_atom(atom.index, atom.resnr, atom.resnr, atom.resname, atomtype);
      output.push_back(new_atom);
      break;
    }
  }
  return output;
}

std::vector<PosResInst> boxtools::makeRestraints(const Box& box, const boxtools::ParameterTable<PosResType>& t){
  std::vector<PosResInst> output;
  if(t.entries.size() == 0) return output;
  for(const auto& atom : box.atoms){
    for(auto& restype : t.entries){
      if(atom.name != restype.name) continue;
      PosResInst new_res(atom.index, restype);
      output.push_back(new_res);
      break;
    }
  }
  return output; 
}

std::vector<BondInst> boxtools::makeBonds(const Box& box, const boxtools::ParameterTable<BondType>& btable)
{
  std::vector<BondInst> bonds;
  if(btable.entries.size() == 0) return bonds;
  int lastRes = -1, natoms = box.atoms.size();
  for(int i = 0; i < natoms; i++){
    for(int j = i+1; j < natoms; j++){
      if(box.atoms[i].resnr != box.atoms[j].resnr) break;
      int pos = isBondInTable(box.atoms[i], box.atoms[j], btable);
      if(pos == -1) continue;
      BondInst new_bond(box.atoms[i].index, box.atoms[j].index, btable.entries[pos]);
      bonds.push_back(new_bond);
    }
  }
  return bonds;
}


std::vector<BondInst> boxtools::makePeriodicBonds(Box& box, const boxtools::ParameterTable<PBCBondType>& pbc_bonds)
{
	std::vector<BondInst> bonds;
  if(pbc_bonds.entries.size() == 0) return bonds;
  Vec3<double> box_dims = {box.boxvec[0][0], box.boxvec[1][1], box.boxvec[2][2]};
  for(const PBCBondType& pbcbond : pbc_bonds.entries){
    double rmax2 = pbcbond.params[0];
    rmax2 *= rmax2;
    auto pbc = pbcbond.pbc;
    auto n1 = pbcbond.i, n2 = pbcbond.j;
    for(std::size_t i = 0; i < box.atoms.size(); i++)
    {
      const Atom& a1 = box.atoms[i];
      for(std::size_t j = i+1; j < box.atoms.size(); j++)
      {
        const Atom& a2 = box.atoms[j];
        if( !((a1.name == n1 && a2.name == n2) || (a1.name == n2 && a2.name == n1)) ) continue;
        double r2 = getPeriodicDistance2(a1, a2, box_dims, pbc);
        if( r2 < rmax2 ){
          BondInst binst;
          binst.i = box.atoms[i].index;
          binst.j = box.atoms[j].index;
          binst.params =  pbcbond.params;
          binst.params[0] = std::sqrt(r2);
          binst.function_type = pbcbond.function_type;
          bonds.push_back(binst);
        }
      }
    }
  }
	return bonds;
}
//returns 1 if one of an atom isn't in the table
bool setAtomParams(Box& box, const boxtools::ParameterTable<AtomType>& atable){
  bool retval = 0;
  for(auto& atom : box.atoms){
    int pos = isAtomInTable(atom, atable);
    if(pos == -1){
      retval = 1;
      continue;
    }
    atom.charge = atable.entries[pos].charge;
    atom.mass = atable.entries[pos].mass;
    atom.type = atable.entries[pos].type;
  }
  return retval;
}
std::vector<AngleInst> boxtools::makeAngles(const Box& box, const boxtools::ParameterTable<AngleType>& atable)
{
  std::vector<AngleInst> angles;
  if( atable.entries.size() == 0) return angles;
  int natoms = box.atoms.size();
  for(int i = 0; i < natoms; i++){
    for(int j = i+1; j < natoms; j++){
      if(box.atoms[i].resnr != box.atoms[j].resnr) break;
      for(int k = j+1; k < natoms; k++){
        if(box.atoms[k].resnr != box.atoms[j].resnr) break;
        int pos;
        int swap = isAngleInTable(box.atoms[i], box.atoms[j], box.atoms[k], atable, pos);
        if(pos == -1) continue;
        if(swap == 0) angles.push_back(AngleInst(box.atoms[i].index, box.atoms[j].index, box.atoms[k].index, atable.entries[pos]));
        else if(swap == 1) angles.push_back(AngleInst(box.atoms[i].index, box.atoms[k].index, box.atoms[j].index, atable.entries[pos]));
        else if(swap == 2) angles.push_back(AngleInst(box.atoms[j].index, box.atoms[i].index, box.atoms[k].index, atable.entries[pos]));
      }
    }
  }
  return angles;
}

std::vector<AngleInst> boxtools::makeAnglesUsingBonds(const Box& box, const boxtools::ParameterTable<AngleType>& atable, const std::vector<BondInst>& bonds){
  std::vector<AngleInst> angles;
  if( atable.entries.size() == 0) return angles;
  int natoms = box.atoms.size();
  for(int i = 0; i < bonds.size(); i++){
    for(int j = i+1; j < bonds.size(); j++){
      auto& bond1 = bonds[i];
      auto& bond2 = bonds[j];
      int a,b,c;
      if(bond1.i == bond2.i){
        a = bond1.j;
        b = bond1.i;
        c = bond2.j;
      }
      else if(bond1.i == bond2.j){
        a = bond1.j;
        b = bond1.i;
        c = bond2.i;
      }
      else if(bond1.j == bond2.i){
        a = bond1.i;
        b = bond1.j;
        c = bond2.j;
      }
      else if(bond1.j == bond2.j){
        a = bond1.i;
        b = bond1.j;
        c = bond2.i;
      }
      else continue; //keep going if there are no shared indices
      //kludge, need more reliable way to go from gro file index to true index
      const Atom& a1 = box.atoms[a-1], a2 = box.atoms[b-1], a3 = box.atoms[c-1];
      int pos;
      int swap = isAngleInTable(a1, a2, a3, atable, pos);
      if(pos == -1) continue;
      if(swap == 0) angles.push_back(AngleInst(a, b, c, atable.entries[pos])); //no index-swapping here since the bonding information is important
    }
  }
  return angles;  
}

std::vector<ExclusionInst> boxtools::makeExclusions(Box& box, const boxtools::ParameterTable<ExclusionType>& t){
  std::vector<ExclusionInst> bonds;
  if(t.entries.size() == 0) return bonds;
  int natoms = box.atoms.size();
  for(int i = 0; i < natoms; i++){
    std::vector<int> indexes;
    for(int j = i+1; j < natoms; j++){
      if(box.atoms[i].resnr != box.atoms[j].resnr) break;
      if(!isExclusionInTable(box.atoms[i], box.atoms[j], t)) continue;
      indexes.push_back(box.atoms[j].index);
    }
    if(indexes.size() > 0){
      ExclusionInst newExcl(box.atoms[i].index, indexes);
      bonds.push_back(newExcl);
    }
  }
  return bonds;  
}

std::vector<ConstraintInst> boxtools::makeConstraints(Box& box, const boxtools::ParameterTable<ConstraintType>& t){
  std::vector<ConstraintInst> constraints;
  if(t.entries.size() == 0) return constraints;
  int natoms = box.atoms.size();
  for(int i = 0; i < natoms; i++){
    for(int j = i+1; j < natoms; j++){
      if(box.atoms[i].resnr != box.atoms[j].resnr) break;
      int pos = isConstraintInTable(box.atoms[i], box.atoms[j], t);
      if(pos == -1) continue;
      auto newconst = ConstraintInst(box.atoms[i].index, box.atoms[j].index, t.entries[pos]);
      constraints.push_back(newconst);
    }
  }
  return constraints;
}


std::vector<int> getIndexesWithResnr(const Box& box, int resnr, int refIdx){
  std::vector<int> ret;
  for(int i = 0; i < box.atoms.size(); i++){
    if(i == refIdx) continue;
    if(box.atoms[i].resnr == resnr) ret.push_back(i);
  }
  return ret;
}

std::vector<Vsite3Inst>  boxtools::makeVsites(Box& box, const boxtools::ParameterTable<Vsite3Type>& t){
  std::vector<Vsite3Inst> vsites;
  if(t.entries.size() == 0) return vsites;
  int lastRes = -1, natoms = box.atoms.size();
  std::vector<int> resIdx;
  for(int i = 0; i < natoms; i++){
    lastRes = box.atoms[i].resnr;
    resIdx = getIndexesWithResnr(box, lastRes, i);
    //std::cin.get(); 
    int nresatoms = resIdx.size();
    for(int j = 0; j < nresatoms; j++){
      for(int k = 0; k < nresatoms; k++){
        if(j == k) continue;
        for(int l = 0; l < nresatoms; l++){
          if(j == l || k ==  l) continue;
          int idx = i;
          int jdx = resIdx[j];
          int kdx = resIdx[k];
          int ldx = resIdx[l];
          auto pos = isVsiteInTable(box.atoms[idx].name, box.atoms[jdx].name, box.atoms[kdx].name, box.atoms[ldx].name, t);
          if(pos == -1) continue;
          auto newvsite = Vsite3Inst(box.atoms[idx].index, box.atoms[jdx].index, box.atoms[kdx].index, box.atoms[ldx].index, t.entries[pos]);
          vsites.push_back(newvsite);
        }
      }
    }
  }
  return vsites;  
} 
