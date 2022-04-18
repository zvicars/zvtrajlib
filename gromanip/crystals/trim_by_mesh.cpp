#include "../actions.hpp"
#include "crystals.hpp"
#include "../../tools/meshtools.hpp"
bool isInsideMesh(std::array<double,3> pos, const Mesh& mesh){
  //needs to be fed a reasonably water-tight mesh to work
  //also naive, slow algorithm so should stick to small meshes
  int count = 0;
  for(auto triangle : mesh.triangles){
    std::array<std::array<double,3>, 3> triangle_pos;
    std::array<double,3> v1, v2, v3;
    for(int i = 0; i < 3; i++){
      auto& vtx = mesh.vertices[triangle.indices[i]];
      triangle_pos[i] = vtx;
    }
    std::array<double,3> inter;
    bool collides = MollerTrombore(pos, {.111, .333, .555}, triangle_pos, inter);
    if(collides) count++;
  }
  if(count%2 == 0) return 0;
  return 1;
}
void boxtools::actions::trimbymesh(GroManipData& data, const std::vector<std::string>& args){
  //argument is a list of volume names
  FANCY_ASSERT(args.size() == 4, "Invalid call to boxtools::actions::trimbymesh(), needs input box, meshfile, mode, and output box name");
  std::string in = args[0], out = args[3];
  std::string meshfile = args[1];
  int mode = std::stoi(args[2]);
  FANCY_ASSERT(mode == 0 || mode == 1, "Invalid mode. Mode 0 is remove all atoms outside of mesh, Mode 1 removes all atoms inside of mesh");
  std::string ag = args[2];
  std::vector<std::string>  new_args = args;
  Box *b1 = data.findBox(in);
  FANCY_ASSERT(b1 != 0, "boxtools::actions::trimbymesh() box not found");
  Box box_out = *b1;
  Box* b_out = data.findBox(out);
  if(b_out == 0){
    b_out = new Box;
  }

  Mesh mesh = getMesh(meshfile);

  std::set<int> resnums;
  for(auto& atom : box_out.atoms){
    auto pos = atom.x;
    if(mode == 0){
      if(isInsideMesh(pos, mesh)) resnums.insert(atom.resnr);
    }
    else if (mode == 1){
      if(!isInsideMesh(pos, mesh)) resnums.insert(atom.resnr);
    }
  }
  std::vector<int> reslist;
  reslist.insert(reslist.begin(), resnums.begin(), resnums.end());
  removeResNumbers(box_out, reslist);
  renumberBox(box_out);
  *b_out = box_out;
  data.addBox(out, b_out); 
  return;
}