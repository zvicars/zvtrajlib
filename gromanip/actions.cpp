#include "actions.hpp"
#include "../tools/pbcfunctions.hpp"
#include "../tools/StringTools.hpp"
#include "../tools/stlmath.hpp"
#include "../tools/cellgrid.hpp"
#include "./add_particle_extra/safe_particle_insert.hpp"
#include <set>
#include <random>

void boxtools::actions::makebox(GroManipData& data, const std::vector<std::string>& args){
    //takes a filename and an output box name
  FANCY_ASSERT(args.size() == 1, "Invalid call to boxtools::actions::makebox(), requires box_name");
  Box* b_out = data.findBox(args[0]);
  if(b_out == 0){
    b_out = new Box;
  }
  data.addBox(args[0], b_out);
  return;
}

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
void boxtools::actions::trimvolumebyresnameinvperiodic(GroManipData& data, const std::vector<std::string>& args){
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
  boxtools::removeResNumbers(box_out, boxtools::getResnrNotWithinVolumePeriodic(box_out, *v1));
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

void boxtools::actions::trimvolumebyatomnameinvperiodic(GroManipData& data, const std::vector<std::string>& args){
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
  boxtools::removeResNumbers(box_out, boxtools::getResnrWithinVolumebyAtomNamePeriodic(box_out, *v1, args[2]));
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
void boxtools::actions::invrotate(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, euler angles (degrees), and an output box name
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::invrotate(), \
  requires input_box_name, euler phi, euler theta, euler psi, output_box_name");
  Vec3<double> angles; 
  for(int i = 0; i < 3; i++){
    angles[i] = std::stod(args[i+1]);
  }
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::invrotate()"); 
  Box* b_out = data.findBox(args[4]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box box_out = *b1;
  boxtools::invrotateEulerAngles(box_out, angles);
  *b_out = box_out;
  data.addBox(args[4], b_out); 
  return;
}
void boxtools::actions::rotate_vector(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, euler angles (degrees), and an output box name
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::rotate_vector(), \
  requires input_box_name, [v1,v2,v3], [w1,w2,w3], output_box_name");
  std::string in = args[0], out = args[3];
  Vec3<double> angles;
  auto v1 = StringTools::stringToVector<double>(args[1]);
  auto v2 = StringTools::stringToVector<double>(args[2]);
  auto a1 = vec2Array<double,3>(v1);
  auto a2 = vec2Array<double,3>(v2);
  a1 = a1 * 1.0/norm2(a1);
  a2 = a2 * 1.0/norm2(a2);
  Box *b1 = data.findBox(in);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::rotate_vector()"); 
  Box* b_out = data.findBox(out);
  if(b_out == 0){
    b_out = new Box;
  }
  Box box_out = *b1;
  //rotation matrix that rotates v1 to v2 
  boxtools::rotateVectorCOM(box_out, a1, a2);
  *b_out = box_out;
  data.addBox(out, b_out); 
  return;
}
void boxtools::actions::flip(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name an axis, and an output box name
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::flip(), requires input_box_name, axis_name, and output_box_name");
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::rotate()"); 
  Box* b_out = data.findBox(args[2]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box box_out = *b1;
  int axis = std::stoi(args[1]);
  boxtools::flipAtoms(box_out, axis);
  *b_out = box_out;
  data.addBox(args[2], b_out); 
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

void boxtools::actions::scale(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a translation vector, and an output box name
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::scale(), requires input_box_name, \
  x factor, y factor, z factor, and output_box_name");
  Vec3<double> pos;
  for(int i = 0; i < 3; i++){
    pos[i] = std::stod(args[i+1]);
  }
  Box *b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::scale()"); 
  Box* b_out = data.findBox(args[4]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box box_out = *b1;
  boxtools::scaleAtoms(box_out, pos);
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
  renumberBox(box_out);
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

void boxtools::actions::deleterandom(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::deleterandom(), needs input box, resname, count, and output box name");
  std::vector<std::string>  new_args = args;
  std::string input_name = args[0], output_name = args[3], resname = args[1];
  int count = std::stoi(args[2]);
  Box *b1 = data.findBox(input_name);
  FANCY_ASSERT(b1 != 0, "boxtools::actions::deleterandom() box not found");
  Box box_out = *b1;
  Box* b_out = data.findBox(output_name);
  if(b_out == 0){
    b_out = new Box;
  }
  int nres = countRes(box_out, resname);
  int deleteCounter = nres - count;
  std::vector<int> vec_resnrs;
  std::set<int> set_resnrs;
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(1, nres);
  while(set_resnrs.size() < deleteCounter){
    int idx = distribution(generator);
    set_resnrs.insert(idx);
  }
  vec_resnrs.insert(vec_resnrs.end(), set_resnrs.begin(), set_resnrs.end());
  removeResNumbers(box_out,vec_resnrs);
  *b_out = box_out;
  data.addBox(output_name, b_out); 
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
  std::cout << exclusiontypes << std::endl;
  std::string constrainttypes = getGMXTaggedParams(filetext, "constrainttypes");
  std::string angletypes = getGMXTaggedParams(filetext,"angletypes");
  std::string angle2types = getGMXTaggedParams(filetext,"angletypes2"); //uses bonded information
  std::string angle3types = getGMXTaggedParams(filetext,"angletypes3"); //uses bonded information  
  std::string dihedraltypes = getGMXTaggedParams(filetext,"dihedraltypes");
  std::string vsitetypes = getGMXTaggedParams(filetext,"vsite3types");

  auto AtomTable = generateTable<BTAtomType>(atomtypes);
  auto PosResTable = generateTable<PosResType>(posrestypes);
  auto BondTable = generateTable<BondType>(bondtypes);
  auto PBCBondTable = generateTable<PBCBondType>(pbcbondtypes);
  auto AngleTable = generateTable<AngleType>(angletypes);
  auto AngleTable2 = generateTable<AngleType>(angle2types);
  auto AngleTable3 = generateTable<AngleType>(angle3types);
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
  std::vector<AngleInst> angles3 = makePBCAnglesUsingBonds(*b1, AngleTable3, bonds);
  angles.insert(angles.end(), angles2.begin(), angles2.end());
  angles.insert(angles.end(), angles3.begin(), angles3.end());
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
//good for adding a single atom to the dataset to balance the charge
void boxtools::actions::addatom(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::addatom(), needs input box, name, resname, position [x,y,z], output box name");
  std::string in = args[0], out = args[4];
  auto pos = StringTools::stringToVector<double>(args[3]);
  std::string atom_name = args[1], resname = args[2];
  Box *b1 = data.findBox(in);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::addatom()");
  Box box_out = *b1;
  Box* b_out = data.findBox(out);
  if(b_out == 0){
    b_out = new Box;
  }
  Atom new_atom;
  new_atom.x = vec2Array<double,3>(pos);
  new_atom.name = atom_name;
  new_atom.resname = resname;
  new_atom.index = box_out.atoms.back().index + 1;
  new_atom.resnr = box_out.atoms.back().resnr + 1;
  new_atom.type = atom_name;
  box_out.atoms.push_back(new_atom);
  *b_out = box_out;
  data.addBox(out, b_out); 
  return;
}

//adds a set of single-atom residues within a given volume
//argument list: input box, volume_name, minimum distance, resname, atomname, quantity, output box
//resname, atomname, and quantity can be repeated multiple times, so (nargs-4)%3 must be 0 with string, string, int args
void boxtools::actions::addparticles(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  int nargs = args.size();
  FANCY_ASSERT(nargs>=7, "Invalid call to boxtools::actions::particles(), need at least 7 arguments");
  FANCY_ASSERT((nargs-4)%3 == 0, "Invalid call to boxtools::actions::particles(), look at code for function arguments");
  std::string in = args[0], out = args[nargs-1];
  std::string volume_name = args[1];
  Volume* v1 = data.findVolume(volume_name);
  FANCY_ASSERT(v1 != 0,  "Failed to find input volume in boxtools::actions::trimvolume()");
  double min_dist = std::stod(args[2]);

  std::vector<std::string> atom_args;
  for(int i = 3; i < nargs-1; i++){
    atom_args.push_back(args[i]);
  }
  int num_atoms = atom_args.size()/3;
  std::vector<std::string> atom_names(num_atoms);
  std::vector<std::string> res_names(num_atoms);
  std::vector<int> atom_counts(num_atoms);
  for(int i=0; i<num_atoms; i++){
    res_names[i] = atom_args[3*i];
    atom_names[i] = atom_args[3*i+1];
    atom_counts[i] = std::stoi(atom_args[3*i+2]);
  }
  std::string atom_name = args[1], resname = args[2];
  Box *b1 = data.findBox(in);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::addatom()");
  Box box_out = *b1;
  Box* b_out = data.findBox(out);
  if(b_out == 0){
    b_out = new Box;
  }
  //meat goes here
  safeInsert s(0.1,v1);
  for(int i = 0; i < num_atoms; i++){
    Atom atom_step;
    atom_step.name = atom_names[i];
    atom_step.resname = res_names[i];
    for(int j = 0; j < atom_counts[i]; j++){
      bool fail = s.getSafePosition(atom_step.x, min_dist);
      FANCY_ASSERT(!fail, "unable to get a safe position, particle density is too high!");
      if(box_out.atoms.size() == 0){
        atom_step.index = 1;
        atom_step.resnr = 1;
      }
      else{
        atom_step.index = box_out.atoms.back().index + 1;
        atom_step.resnr = box_out.atoms.back().resnr + 1;
      }
      box_out.atoms.push_back(atom_step);
    }
  }
  box_out.hasNamedAtoms=1;
  *b_out = box_out;
  data.addBox(out, b_out); 
  return;
}

void boxtools::actions::printindicesnear(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 5, 
  "Invalid call to boxtools::actions::printindicesnear(), needs in, [n1,n2,n3], [n1,n2,n3], distance, outfile");
  std::string in = args.front();
  double distance_thresh = std::stod(args[3]);
  std::string filename = args[4];
  std::string n1 = args[1], n2 = args[2];
  std::vector<std::string> names1, names2;
  names1 = StringTools::stringToVector<std::string>(n1);
  names2 = StringTools::stringToVector<std::string>(n2);
  Box *b1 = data.findBox(in);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::setboxsize()");
  Box box = *b1;
  std::set<int> idx_out;
  std::vector<int> idx1, idx2;
  for(int i = 0; i < box.atoms.size(); i++){
    for(auto name : names1){
      if(box.atoms[i].name == name){
        idx1.push_back(i);
        break;
      }
    }
    for(auto name : names2){
      if(box.atoms[i].name == name){
        idx2.push_back(i);
        break;
      }
    }
  }
  for(auto i : idx1){
    for(auto j : idx2){
      auto pos1 = box.atoms[i].x;
      auto pos2 = box.atoms[j].x;
      double distance = norm2(pos2-pos1);
      if(distance < distance_thresh){
        idx_out.insert(box.atoms[i].index);
      }
    }
  }
  std::ofstream ofile(filename);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for boxtools::actions::printindicesnear()");
  for(auto idx : idx_out){
    ofile << idx << "\n";
  }
  ofile.close();
  return;
}

void boxtools::actions::hollow(GroManipData& data, const std::vector<std::string>& args){
  //this can break up residues, so use judiciously
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::hollow(), \
  requires input_box_name, atom_name, num_neighbors, cutoff distance, and output_box_name");
  Box* b_out = data.findBox(args[3]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::hollow()");
  Box box_out = *b1;
  double cutoff_distance = std::stod(args[2]);
  int num_neighbors = std::stoi(args[1]);
  CellGrid c1;
  auto box_vec = box_out.boxvec;
  Vec3<double> box_vector;
  double box_volume = 0.0;
  for(int i=0; i<3; i++){
    box_vector[i] = box_vec[i][i];
    box_volume *= box_vec[i][i];
  }
  if(box_volume == 0.0){
    for(int i=0; i<3; i++) box_vector[i] = 100; //using large but reasonable number to ensure pbc's don't come into play
  }
  std::set<int> index_list;
  c1.reset(cutoff_distance, box_vector);
  for(int i = 0; i < box_out.atoms.size(); i++){
    c1.addIndexToGrid(i, box_out.atoms[i].x);
  }
  for(int i = 0; i < box_out.atoms.size(); i++){
    auto pos1 = box_out.atoms[i].x;
    auto indices = c1.getNearbyIndices(box_out.atoms[i].x);
    int counter=0;
    for(auto index : indices){
      auto pos2 = box_out.atoms[index].x;
      if(getDistance(pos2, pos1, box_vector) < cutoff_distance){
        counter++;
      }
    }
    counter--; //accounting for self-interaction
    if(counter > num_neighbors) index_list.insert(box_out.atoms[i].resnr);
  }
  std::vector<int> index_list2; 
  index_list2.insert(index_list2.begin(), index_list.begin(), index_list.end());
  removeResNumbers(box_out, index_list2);
  *b_out = box_out;
  data.addBox(args[3], b_out);
  return; 
}

void boxtools::actions::wrap(GroManipData& data, const std::vector<std::string>& args){
  //this can break up residues, so use judiciously
  FANCY_ASSERT(args.size() == 2, "Invalid call to boxtools::actions::wrap(), \
  requires input_box_name, output_box_name");
  Box* b_out = data.findBox(args[1]);
  if(b_out == 0){
    b_out = new Box;
  }
  Box *b1;
  b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::wrap()");
  Box box_out = *b1;
  wrapPBC(box_out);
  *b_out = box_out;
  data.addBox(args[1], b_out);
  return; 
}