#pragma once
#include "../MarchingCubesInterface.hpp"
#include <string>
#include <vector>
#include <array>
#include <map>
class ObjFile{
public:
  int load(std::string file_text);
  int save(std::string& file_text);
  Mesh getMesh();
protected:
  void createPerVertexNormals();
  void readVertex(std::string line);
  void readVertexNormal(std::string line);
  void readFace(std::string line);
private:
  std::vector<std::array<double,3> > vertices_;
  std::vector<int> vertex_normals_;
  std::vector<std::array<int,3> > faces_, face_vtx_norms_;
  std::vector<std::array<double, 3> > vertex_normals_table_;
};