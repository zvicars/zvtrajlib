#pragma once
#include <array>
#include <vector>

template <class T>
using Vec3 = std::array<T, 3>;

struct Triangle{
    Vec3<int> indices; //contains either the real positions of the vertices or the indices, depending on whether you use an int or a double
    Vec3<double> normal; //contains the normal vector of a triangle
};
struct Mesh {
    std::vector<Vec3<double> > vertices;
    std::vector<Vec3<double> > normals;
    std::vector<Triangle> triangles;
    int nvtx, ntri;
};

struct MeshFile{
public:
  virtual int load(std::string file_text)=0;
  virtual int save( const Mesh& mesh, std::string& file_text)=0;
  virtual Mesh getMesh(){
      return meshInternal_;
  }
protected:
  Mesh meshInternal_;
};