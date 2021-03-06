#pragma once
#include "../MarchingCubesInterface.hpp"
#include <string>
#include <vector>
#include <array>
class StlFile{
public:
  virtual int load(std::string file_text);
  virtual int save( const Mesh& mesh, std::string& file_text);
  virtual Mesh getMesh();
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