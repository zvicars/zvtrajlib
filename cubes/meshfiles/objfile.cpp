#include "objfile.hpp"
#include <sstream>

void ObjFile::readVertex(std::string line){
  std::stringstream ss(line);
  std::string tag;
  std::array<double, 3> entry;
  ss >> tag >> entry[0] >> entry[1] >> entry[2];
  if(ss.fail()){
    std::cout << "Bad data in obj file." << std::endl;
    return;
  }
  vertices_.push_back(entry);
  return;
}
void ObjFile::readVertexNormal(std::string line){
  std::stringstream ss(line);
  std::string tag;
  std::array<double, 3> entry;
  ss >> tag >> entry[0] >> entry[1] >> entry[2];
  std::cout << tag << " " << entry[0] << " " << entry[1] << " " << entry[2] << std::endl;
  if(ss.fail()){
    std::cout << "Bad data in obj file." << std::endl;
    return;
  }
  vertex_normals_table_.push_back(entry);
  return;
}
void ObjFile::readFace(std::string line){
  std::stringstream ss(line);
  std::string tag;
  std::array<std::string, 3> entry;
  ss >> tag >> entry[0] >> entry[1] >> entry[2];
  if(ss.fail()){
    std::cout << "Bad data in obj file." << std::endl;
    return;
  }
  std::array<int, 3> face, normal;
  for(int i = 0; i < 3; i++){
    std::string vs, ts, ns;
    int vi, ti, ni;
    int pos1 = entry[i].find('/');
    int pos2 = entry[i].rfind('/');
    if(pos1 == std::string::npos || pos1 == pos2){
      std::cout << "Bad face in obj file" << std::endl;
      return;
    }
    vs = entry[i].substr(0, pos1);
    ts = entry[i].substr(pos1+1, pos2-pos1-1);
    ns = entry[i].substr(pos2+1);
    if(vs.size() == 0){
      std::cout << "Bad face in obj file." << std::endl;
      return;
    }
    if(ts.size() == 0){
      ts = -1;
    }
    if(ns.size() == 0){
      std::cout << "No normal found for face." << std::endl;
      return;      
    }
    face[i] = std::stoi(vs)-1;
    normal[i] = std::stoi(ns)-1;
  }
  std::cout << face[0] << " " << face[1] << " " << face[2] << "\n";
  std::cout << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
  faces_.push_back(face);
  face_vtx_norms_.push_back(normal);
  return;
}
int ObjFile::load(std::string file_text){
  std::stringstream ss(file_text);
  std::string line;
  bool objectFound = 0;
  while(std::getline(ss, line)){
    std::string tag;
    std::stringstream ss2(line);
    ss2 >> tag;
    if(tag == "o"){
      if(objectFound == 1) return 2;
      objectFound == 1;
    }
    else if(tag == "v"){
      readVertex(line);
    }
    else if(tag == "vn"){
      readVertexNormal(line);  
    }
    else if(tag == "f"){
      readFace(line);
    }
    else{
      continue;
    }
  }

  createPerVertexNormals();
  return 0;
}
int ObjFile::save(std::string& file_text){
  return 1;
}

//fills a vector that is the same size as the vertices with the vertex information contained in the faces section
void ObjFile::createPerVertexNormals(){
  vertex_normals_.resize(vertices_.size(), -1);
  for(int i = 0; i < faces_.size(); i++){
    auto face = faces_[i];
    auto normal = face_vtx_norms_[i];
    for(int j = 0; j < 3; j++){
      if(vertex_normals_[face[j]] == -1){
        vertex_normals_[face[j]] = normal[j];
      }
      else if(vertex_normals_[face[j]] != normal[j]){
        std::cout << "Vertex normals disagree with one another. Please recheck your obj file to ensure proper output.\n";
      }
    }
  }

  for(int i = 0; i < vertex_normals_.size(); i++){
    if(vertex_normals_[i] == -1){
      std::cout << "Not all vertices have normal vectors. Please recheck file." << std::endl;
      vertex_normals_[i] = 0; 
    }
  }

  return;
}
Mesh ObjFile::getMesh(){
  Mesh m1;
  std::vector<std::array<double,3> > real_normals(vertices_.size());
  std::vector<Triangle> triangles;
  for(int i = 0; i < vertex_normals_.size(); i++){
    real_normals[i] = vertex_normals_table_[vertex_normals_[i]];
    std::cout << "Real normal line" << std::endl;
    std::cout << real_normals[i][0] << " " << real_normals[i][1] << " " << real_normals[i][2] << "\n";
  }
  for(int i = 0; i < faces_.size(); i++){
    Triangle tri;
    tri.indices = faces_[i];
    tri.normal = vertex_normals_table_[face_vtx_norms_[i][0]];
    triangles.push_back(tri);
  }

  m1.vertices = vertices_;
  m1.normals  = real_normals;
  m1.triangles = triangles;
  m1.ntri = triangles.size();
  m1.nvtx = vertices_.size();

  return m1;
}
