#include "actions.hpp"
#include "../tools/pbcfunctions.hpp"
void boxtools::actions::merge(GroManipData& data, const std::vector<std::string>& args){
  //takes two box names and an output box name
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::merge(), requires input_box_name, input_box_name, and output_box_name");
  Box* b_out = data.findBox(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1, *b2;
  b1 = data.findBox(args[0]);
  b2 = data.findBox(args[1]);
  FANCY_ASSERT(b1 != 0 && b2 != 0, "Failed to find one of the specified boxes in boxtools::actions::merge()");
  Box box_out = boxtools::mergeBox(*b1, *b2);
  *b_out = box_out;
  data.addBox(args[2], b_out);
  return;
}

void boxtools::actions::eraseres(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a resname, and the output box
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::eraseres(), requires input_box_name, resname, and output_box_name");
  Box* b_out = data.findBox(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::eraseres()");
  Box box_out = *b1;
  boxtools::deleteRestype(box_out, args[1]);
  *b_out = box_out;
  data.addBox(args[2], b_out);
  return;
}

void boxtools::actions::keepres(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a resname, and the output box
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::keepres(), requires input_box_name, resname, and output_box_name");
  Box* b_out = data.findBox(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::eraseres()");
  Box box_out = *b1;
  boxtools::deleteNotRestype(box_out, args[1]);
  *b_out = box_out;
  data.addBox(args[2], b_out);
  return;
}

void boxtools::actions::duplicate(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name and an output name
  FANCY_ASSERT(args.size() == 2, "Invalid call to boxtools::actions::duplicate(), requires input_box_name and output_box_name");
  Box* b_out = data.findBox(args[1]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::eraseres()");
  Box box_out = *b1;
  *b_out = box_out;
  data.addBox(args[2], b_out);
  return;
}
void boxtools::actions::trimvolumebyresname(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a volume name, and an output name
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::trimvolume(), requires input_box_name, volume_name and output_box_name");
  Box* b_out = data.findBox(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  Volume *v1;
  b1 = data.findBox(args[0]);
  v1 = data.findVolume(args[1]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::trimvolume()");
  FANCY_ASSERT(v1 != 0,  "Failed to find input volume in boxtools::actions::trimvolume()");
  Box box_out = *b1;
  boxtools::removeResNumbers(box_out, boxtools::getResnrWithinVolume(box_out, *v1));
  *b_out = box_out;
  data.addBox(args[2], b_out);
  return; 
}

void boxtools::actions::trimvolumebyresnameinv(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a volume name, and an output name
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::trimvolume(), requires input_box_name, volume_name and output_box_name");
  Box* b_out = data.findBox(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  Volume *v1;
  b1 = data.findBox(args[0]);
  v1 = data.findVolume(args[1]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::trimvolume()");
  FANCY_ASSERT(v1 != 0,  "Failed to find input volume in boxtools::actions::trimvolume()");
  Box box_out = *b1;
  boxtools::removeResNumbers(box_out, boxtools::getResnrNotWithinVolume(box_out, *v1));
  *b_out = box_out;
  data.addBox(args[2], b_out);
  return; 
}

void boxtools::actions::trimvolumebyatomname(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a volume name, and an output name
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::trimvolumebyatomname(), \
  requires input_box_name, volume_name, atom_name,  and output_box_name");
  Box* b_out = data.findBox(args[3]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  Volume *v1;
  b1 = data.findBox(args[0]);
  v1 = data.findVolume(args[1]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::trimvolumebyatomname()");
  FANCY_ASSERT(v1 != 0,  "Failed to find input volume in boxtools::actions::trimvolumebyatomname()");
  Box box_out = *b1;
  boxtools::removeResNumbers(box_out, boxtools::getResnrNotWithinVolumebyAtomName(box_out, *v1, args[2]));
  *b_out = box_out;
  data.addBox(args[3], b_out);
  return; 
}
void boxtools::actions::trimvolumebyatomnameinv(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a volume name, and an output name
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::trimvolumebyatomname(), \
  requires input_box_name, volume_name, atom_name,  and output_box_name");
  Box* b_out = data.findBox(args[3]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  Volume *v1;
  b1 = data.findBox(args[0]);
  v1 = data.findVolume(args[1]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::trimvolumebyatomname()");
  FANCY_ASSERT(v1 != 0,  "Failed to find input volume in boxtools::actions::trimvolumebyatomname()");
  Box box_out = *b1;
  boxtools::removeResNumbers(box_out, boxtools::getResnrWithinVolumebyAtomName(box_out, *v1, args[2]));
  *b_out = box_out;
  data.addBox(args[3], b_out);
  return; 
}
void boxtools::actions::rotate(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, euler angles (degrees), and an output box name
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::rotate(), \
  requires input_box_name, euler phi, euler theta, euler psi, output_box_name");
  Vec3<double> angles; 
  for(int i = 0; i < 3; i++){
    angles[i] = std::stod(args[i+1]);
  }
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::rotate()"); 
  Box* b_out = data.findBox(args[4]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box box_out = *b1;
  boxtools::rotateEulerAngles(box_out, angles);
  *b_out = box_out;
  data.addBox(args[4], b_out); 
  return;
}
void boxtools::actions::flip(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name an axis, and an output box name
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::flip(), requires input_box_name, axis_name, and output_box_name");
  return;
}
void boxtools::actions::translate(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a translation vector, and an output box name
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::translate(), requires input_box_name, \
  x distance, y distance, z distance, and output_box_name");
  Vec3<double> pos;
  for(int i = 0; i < 3; i++){
    pos[i] = std::stod(args[i+1]);
  }
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::rotate()"); 
  Box* b_out = data.findBox(args[4]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box box_out = *b1;
  boxtools::translateAtoms(box_out, pos);
  *b_out = box_out;
  data.addBox(args[4], b_out); 
  return;
}
void boxtools::actions::loadgro(GroManipData& data, const std::vector<std::string>& args){
  //takes a filename and an output box name
  FANCY_ASSERT(args.size() == 2, "Invalid call to boxtools::actions::loadgro(), requires filename and output_box_name");
  Box* b_out = data.findBox(args[1]);
  if(b_out == 0){
    b_out = new Box;
  }
  readGRO(args[0], *b_out);
  data.addBox(args[1], b_out); 
  return;
}
void boxtools::actions::writegro(GroManipData& data, const std::vector<std::string>& args){
  //takes a an input box name and a filename
  FANCY_ASSERT(args.size() == 2, "Invalid call to boxtools::actions::writegro(), requires input_box_name and filename");
  Box* b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::writegro()"); 
  writeGRO(args[1], b1);
  return;
}
void boxtools::actions::createprimitivevolume(GroManipData& data, const std::vector<std::string>& args){
  //first argument specifies the type of volume, subsequent arguments are fed to the appropriate constructor
  FANCY_ASSERT(args.size() >= 1, "Invalid call to boxtools::actions::loadvolume(), requires volume_type and appropriate arguments for volume");
  FANCY_ASSERT(data.findVolume(args[0]) == 0, "A volume with the specified name already exists.");
  std::vector<std::string> new_args;
  new_args.insert(new_args.begin(), args.begin()+2, args.end());
  auto v1 = createPrimitiveVolume(args[1], new_args);
  data.addVolume(args[0], v1);
  return;
}
void boxtools::actions::createunionvolume(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() >= 1, "Invalid call to boxtools::actions::createunionvolume()");
  std::vector<std::string>  new_args = args;
  new_args.erase(new_args.begin());
  FANCY_ASSERT(data.findVolume(args[0]) == 0, "A volume with the specified name already exists.");
  Volume* v1 = createUnionVolume(data, new_args);
  data.addVolume(args[0], v1);
  return;
}
void boxtools::actions::shrinkwrap(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 2, "Invalid call to boxtools::actions::shrinkwrap(), needs input and output box name");
  std::vector<std::string>  new_args = args;
  new_args.erase(new_args.begin());
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::shrinkwrap()");
  Box box_out = *b1;
  Box* b_out = data.findBox(args[1]);
  if(b_out == 0){
    b_out = new Box;
  }
  shrinkWrap(box_out);
  *b_out = box_out;
  data.addBox(args[1], b_out); 
  return;
}

void boxtools::actions::relabelatom(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::relabelatom(), needs input box, old name, new name, and output box name");
  std::vector<std::string>  new_args = args;
  new_args.erase(new_args.begin());
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "boxtools::actions::relabelatom()");
  Box box_out = *b1;
  Box* b_out = data.findBox(args[3]);
  if(b_out == 0){
    b_out = new Box;
  }
  relabelAtom(box_out, args[1], args[2]);
  *b_out = box_out;
  data.addBox(args[3], b_out); 
  return;
}

void boxtools::actions::relabelres(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::relabelres(), needs input box, old name, new name, and output box name");
  std::vector<std::string>  new_args = args;
  new_args.erase(new_args.begin());
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::relabelres()");
  Box box_out = *b1;
  Box* b_out = data.findBox(args[3]);
  if(b_out == 0){
    b_out = new Box;
  }
  relabelRes(box_out, args[1], args[2]);
  *b_out = box_out;
  data.addBox(args[3], b_out); 
  return;
}
void boxtools::actions::respattern(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::relabelres(), needs input box, number, and output box name");
  std::vector<std::string>  new_args = args;
  new_args.erase(new_args.begin());
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::relabelres()");
  Box box_out = *b1;
  Box* b_out = data.findBox(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  int pattern = std::stoi(args[1]);
  int counter = 0;
  int rescounter = 1;
  for(auto& atom : box_out.atoms){
    if(counter % pattern == 0) rescounter++;
    atom.resnr = rescounter;
    counter++;
  }
  *b_out = box_out;
  data.addBox(args[2], b_out); 
  return;
}
void boxtools::actions::setboxsize(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::setboxsize(), needs input box, x, y, z, output box name");
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::setboxsize()");
  Box box_out = *b1;
  Box* b_out = data.findBox(args[4]);
  if(b_out == 0){
    b_out = new Box;
  }
  Vec3<double> box_vec;
  box_vec.fill(0.0);
  for(int i = 0; i < 3; i++){
    box_vec[i] = std::stod(args[i+1]);
  }
  setBoxSize(box_out, box_vec);
  *b_out = box_out;
  data.addBox(args[4], b_out); 
  return;
}

//assumes entire box is one big molecule
void boxtools::actions::outputmolecule(GroManipData& data, const std::vector<std::string>& args){
  //a molecule definition has the following components
  //atomtype definition : maps atom names to atom types and specifies mass/charge
  //[ molecule ] provides a name for the molecule (mandatory)
  //[ atomtypes ] maps gro file atom names to atom types and mass/charge (mandatory)
  //[ bondtypes ] maps atom name pairs to  bonded potentials, operates on atoms with the same resnr
  //[ pbcbondtypes ] maps atom name pairs to the bonded potentials, works between resnumbers and calculates appropriate bond distance across pbcs
  //[ settlestypes ] specifies settles directive and the name of the oxygen atom to apply it to
  //[ exclusiontypes ] specifies atom name pairs for which to generate exclusions
  //[ angletypes ] maps atom name triplets to angle potentials, operates on atoms with the same resnr
  //[ vsitetypes ] generates vsite3 entries where the first entry is the virtual atom
  std::string filename, ofilename;
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::outputmolecule(), needs input box, ref file name, and output file");
  filename = args[1];
  ofilename = args[2];
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::setboxsize()");

  std::string filetext = loadFileIntoString(filename);
  std::string molname = getGMXTaggedParams(filetext, "molecule");
  std::string atomtypes = getGMXTaggedParams(filetext,"atomtypes");
  std::string posrestypes = getGMXTaggedParams(filetext,"position_restraints");
  std::string bondtypes = getGMXTaggedParams(filetext,"bondtypes");
  std::string pbcbondtypes = getGMXTaggedParams(filetext,"pbcbondtypes");
  std::string exclusiontypes = getGMXTaggedParams(filetext, "exclusiontypes");
  std::string constrainttypes = getGMXTaggedParams(filetext, "constrainttypes");
  std::string angletypes = getGMXTaggedParams(filetext,"angletypes");
  std::string angle2types = getGMXTaggedParams(filetext,"angletypes2"); //uses bonded information
  std::string dihedraltypes = getGMXTaggedParams(filetext,"dihedraltypes");
  std::string vsitetypes = getGMXTaggedParams(filetext,"vsite3types");

  auto AtomTable = generateTable<BTAtomType>(atomtypes);
  auto PosResTable = generateTable<PosResType>(posrestypes);
  auto BondTable = generateTable<BondType>(bondtypes);
  auto PBCBondTable = generateTable<PBCBondType>(pbcbondtypes);
  auto AngleTable = generateTable<AngleType>(angletypes);
  auto AngleTable2 = generateTable<AngleType>(angle2types);
  auto ExclusionTable = generateTable<ExclusionType>(exclusiontypes);
  auto ConstraintTable = generateTable<ConstraintType>(constrainttypes);
  auto VsiteTable = generateTable<Vsite3Type>(vsitetypes);

  std::vector<AtomInst> atoms = makeAtoms(*b1, AtomTable);
  std::vector<PosResInst> restraints = makeRestraints(*b1, PosResTable); 
  std::vector<BondInst> bonds = makeBonds(*b1, BondTable);
  std::vector<BondInst> pbc_bonds = makePeriodicBonds(*b1, PBCBondTable);
  bonds.insert(bonds.end(), pbc_bonds.begin(), pbc_bonds.end());
  std::vector<AngleInst> angles = makeAngles(*b1, AngleTable);
  std::vector<AngleInst> angles2 = makeAnglesUsingBonds(*b1, AngleTable2, bonds);
  angles.insert(angles.end(), angles2.begin(), angles2.end());
  std::vector<ExclusionInst> exclusions = makeExclusions(*b1, ExclusionTable);
  std::vector<ConstraintInst> constraints = makeConstraints(*b1, ConstraintTable);
  std::vector<Vsite3Inst> vsites = makeVsites(*b1, VsiteTable);
  std::ofstream ofile(ofilename);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for molecule bond data");
  
  ofile << "[ moleculetype ]" << "\n";
  ofile << molname << "\n";

  int natoms = atoms.size();
  if(natoms > 0){
    ofile << "[ atoms ]" << "\n";
    for(int i = 0; i < natoms; i++){
      ofile << atoms[i].print() << "\n";
    }
    ofile << std::endl;
  }
  int nrest = restraints.size();
  if(nrest > 0){
    ofile << "[ position_restraints ]" << "\n";
    for(int i = 0; i < nrest; i++){
      ofile << restraints[i].print() << "\n";
    }
    ofile << std::endl;
  }
  int nbonds = bonds.size();
  if(nbonds > 0){
    ofile << "[ bonds ]" << "\n";
    for(int i = 0; i < nbonds; i++){
      ofile << bonds[i].print() << "\n";
    }
  ofile << std::endl;
  }
  int nexcl = exclusions.size();
  if(nexcl > 0){
    ofile << "[ exclusions ]" << "\n";
    for(int i = 0; i < nexcl; i++){
      ofile << exclusions[i].print() << "\n";
    }
    ofile << std::endl;
  }
  int nconst = constraints.size();
  if(nconst > 0){
    ofile << "[ constraints ]" << "\n";
    for(int i = 0; i < nconst; i++){
      ofile << constraints[i].print() << "\n";
    }
    ofile << std::endl;
  }

  int nangles = angles.size();
  if(nangles > 0){
    ofile << "[ angles ]" << "\n";
    for(int i = 0; i < nangles; i++){
      ofile << angles[i].print() << "\n";
    }
    ofile << std::endl;
  }
  int nv = vsites.size();
  if(nv > 0){
    ofile << "[ virtual_sites3 ]" << "\n";
    for(int i = 0; i < nv; i++){
      ofile << vsites[i].print() << "\n";
    }
    ofile << std::endl;
  }
  ofile.close();
  return;

}

void boxtools::actions::deleteoverlapping(GroManipData& data, const std::vector<std::string>& args){
  //takes two box names and an output box name
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::deleteoverlapping(), requires input_box_name, input_box_name, threshold, and output_box_name");
  Box* b_out = data.findBox(args[3]);
  double thresh = std::stod(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1, *b2;
  b1 = data.findBox(args[0]);
  b2 = data.findBox(args[1]);
  FANCY_ASSERT(b1 != 0 && b2 != 0, "Failed to find one of the specified boxes in boxtools::actions::deleteoverlapping()");
  Box box_out = *b1;
  std::vector<int> resnums_to_delete;
  for(auto& atoms : box_out.atoms){
    for(const auto& atoms2 : b2->atoms){
      double dist = 0.0;
      for(int i = 0; i < 3; i++){
        dist += (atoms.x[i] - atoms2.x[i])*(atoms.x[i] - atoms2.x[i]);
      }
      dist = std::sqrt(dist);
      if( dist < thresh) resnums_to_delete.push_back(atoms.resnr);
    }
  }
  removeResNumbers(box_out, resnums_to_delete);
  *b_out = box_out;
  data.addBox(args[3], b_out);
  return;
} 


void boxtools::actions::pbccorrect(GroManipData& data, const std::vector<std::string>& args){
  //takes two box names and an output box name
  FANCY_ASSERT(args.size() == 2, "Invalid call to boxtools::actions::pbccorrect(), requires input_box_name and output_box_name");
  Box* b_out = data.findBox(args[1]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find specified box in boxtools::actions::pbccorrect()");
  Box box_out = *b1;
  std::array<double,3> box_dims;
  for(int i = 0; i < 3; i++) box_dims[i] = box_out.boxvec[i][i];
  for(auto& atom : box_out.atoms){
    placeInsideBox(atom.x, box_dims);
  }
  *b_out = box_out;
  data.addBox(args[1], b_out);
  return;
}