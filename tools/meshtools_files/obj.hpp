#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include "meshfiles.hpp"
class ObjFile : public MeshFile{
public:
  int load(std::string file_text);
  int save( const Mesh& mesh, std::string& file_text);
protected:
  void createPerVertexNormals();
  void readVertex(std::string line);
  void readVertexNormal(std::string line);
  void readFace(std::string line);
  void setMesh();
  std::vector<std::array<double,3> > vertices_;
  std::vector<int> vertex_normals_;
  std::vector<std::array<int,3> > faces_, face_vtx_norms_;
  std::vector<std::array<double, 3> > vertex_normals_table_;
};