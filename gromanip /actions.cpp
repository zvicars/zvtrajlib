#include "actions.hpp"
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
void boxtools::actions::trimvolumebyatomname(GroManipData& data, const std::vector<std::string>& args){
  //takes a box name, a volume name, and an output name
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::trimvolumebyatomname(), \
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
  Box* b_out = data.findBox(args[3]);
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
  FANCY_ASSERT(args.size() == 3, "Invalid call to boxtools::actions::writegro(), requires input_box_name and filename");
  Box* b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::writegro()"); 
  writeGRO(args[1], b1);
  return;
}
void boxtools::actions::createprimitivevolume(GroManipData& data, const std::vector<std::string>& args){
  //first argument specifies the type of volume, subsequent arguments are fed to the appropriate constructor
  FANCY_ASSERT(args.size() >= 1, "Invalid call to boxtools::actions::loadvolume(), requires volume_type and appropriate arguments for volume");
  FANCY_ASSERT(data.findVolume(args[0]) == 0, "A volume with the specified name already exists.");
  std::vector<std::string>  new_args = args;
  new_args.erase(new_args.begin());
  auto v1 = createPrimitiveVolume(args[0], new_args);
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
