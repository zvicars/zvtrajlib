#pragma once
#include <vector>
#include <array>
#include <string>
#include <map>
template <class T>
using Vec3 = std::array<T, 3>;

template <class T>
using Vec = std::vector<T>;

struct Atom{ //many of these properties are optional
  Vec3<double> x, v;
  double mass, charge;
  int resnr, index;
  std::string name, type, resname;
};

struct AtomType{
  std::string name, type;
  double mass, charge;
  int nbfunc;
  Vec<double> nonbonded_params;
};

struct AtomDef{
  int nr, resnr, cgnr, charge, mass; //mass and charge can default to atomtype definition
  std::string type, res, name; //type is the pre-defined atom type, name is the alias that will be in the gro file
  int nbfunc;
  double V, W; //based on gmx definition, will look at defaults for nbfunc
};

struct Bond{
  int i, j;
  std::vector<double> params;
  int function_type;
  double (*calculate)(const Vec3<double>&, const Vec3<double>&, const Vec<double>&) = 0;
};

struct Pair{
  int i, j;
  std::vector<double> params;
  int function_type;
  double (*calculate)(const Vec3<double>&, const Vec3<double>&, const Vec<double>&) = 0;
};

struct Angle{
  int i, j, k;
  std::vector<double> params;
  int function_type;
  double (*calculate)(const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, Vec<double>&) = 0;
};


struct Dihedral{
  int i, j, k, l;
  std::vector<double> params;
  int function_type;
  double (*calculate)(const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, Vec<double>&) = 0;
};

struct MolType{
  //; nr type resnr residue atom cgnr charge mass
  Vec<AtomDef> atoms; //maps to registered set of atom types
  Vec<Bond> bonds;
  Vec<Angle> angles;
  Vec<Dihedral> dihedrals;
  Vec<Pair> pairs; // 1-4 pairs, can be derived from bonds but will most likely be read from top file
};

struct TopologyInfo{
  Vec<AtomType> atomtypes;
  Vec<MolType> molecules;
  Vec<int> counts;
};
struct IndexInfo{
  std::map<std::string, Vec<int> > indexes;
};

class Box{
public:
  Vec<Atom> atoms;
  IndexInfo idxinfo;
  TopologyInfo topinfo;
  Vec3<Vec3<double> > boxvec;
  double time;
  int frame, frame_counter = -1;
  bool hasNamedAtoms = 0, hasTopologyInfo = 0, hasIndexes = 0, hasVelocities = 0;
  Vec3<double> getAtomPosition(int idx) const{
    return atoms[idx].x;
  }
};