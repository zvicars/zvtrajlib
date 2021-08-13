#include "MarchingCubesInterface.hpp"
#include "../extern/eigen-3.3.9/Eigen/Eigen"
#include "../extern/eigen-3.3.9/unsupported/Eigen/NonLinearOptimization"
#include <set>

Eigen::Matrix3d getVertexTransformationMatrix(const Mesh& mesh, int idx){
  auto vertex_in = mesh.vertices[idx];
  auto normal_in = mesh.normals[idx];
  Eigen::Vector3d vertex, normal;
  Eigen::Vector3d unit_z;
  unit_z << 0, 0, 1;
  for(int i= 0; i < 3; i++){
    vertex(i) = vertex_in[i];
    normal(i) = normal_in[i];
  }
  normal.normalize();
  auto v = normal.cross(unit_z);
  double c = normal.dot(unit_z);
  double s = v.norm();
  Eigen::Matrix3d kmat, eye_mat, rotation_matrix;
  kmat << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
  rotation_matrix = Eigen::Matrix3d::Identity() + kmat + Eigen::Matrix3d::Ones()*(kmat.dot(kmat)*((1-c)/(s*s)));
  return rotation_matrix;
}

std::vector<int> getNeighboringIndices(const Mesh& mesh, int idx){
  std::vector<int> indices;
  std::array<int, 3> indices_step;
  for(int i = 0; i < mesh.triangles.size(); i++){
    bool step_flag = 0;
    indices_step = mesh.triangles[i].indices;
    for(int j = 0; j < 3; j++){
      if(indices_step[i] == idx){
        step_flag = 1;
      }
    }
    for(int j = 0; j < 3; j++){
      if(indices_step[i] != idx) indices.push_back(indices_step[i]);
    }
  }
  return indices;
}

inline std::vector<int> getNeighborsWithinRange(const Mesh& mesh, int idx, int neighbors){
  std::vector<int> indices(1, idx);
  std::vector<int> last_indices(1, idx);
  std::vector<int> new_indices;
  for(int i = 0; i < neighbors; i++){
    last_indices = new_indices;
    new_indices.clear();
    for(int j = 0; j < last_indices.size(); j++){
      std::vector<int> new_index_step = getNeighboringIndices(mesh, j);
      new_indices.insert(new_indices.end(), new_index_step.begin(), new_index_step.end());
    }
    for(int j = 0; j < new_indices.size(); j++){
      bool containsIndex = 0;
      for(int k = 0; k < indices.size(); k++){
        if(new_indices[j] == indices[k]) containsIndex = 1;
      }
      if(containsIndex){
        new_indices.erase(new_indices.begin() + j);
        j--;
      }
    }
    indices.insert(indices.end(), new_indices.begin(), new_indices.end());
  }
  return indices;
}

struct Functor2Form{
  Eigen::MatrixXd data;
  Functor2Form(Eigen::MatrixXd data_in){
    data = data_in;
    return;
  }
  int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec){
    assert(b.size() == 3);
    fvec.resize(data.rows()); //make sure the f vector has the same number of rows as the data vector
    for(int i = 0; i < fvec.size(); i++){
      fvec[i] = 0.5*b(0)*data(i,0)*data(i,1) + b(1)*data(i,0)*data(i,1) + 0.5*b(2)*data(i,1)*data(i,1) - data(i,2);
    }
    return 0;
  }
  int df(const Eigen::VectorXd &b, Eigen::MatrixXd &fjac)
  {
    assert(b.size() == 3);
    fjac.resize(data.rows(), 3);
    for(int i = 0; i < data.rows(); i++){
      Eigen::Vector3d jac_row;
      jac_row << 0.5*data(i,0)*data(i,0), data(i,0)*data(i,1), 0.5*data(i,1)*data(i,1);
      fjac.row(i) = jac_row;
    }
    return 0;
  }
};

Eigen::Matrix2d get2ndFormTensor(const Mesh& mesh, int idx){
  Eigen::Matrix2d eval;
  std::vector<int> included_indices = getNeighborsWithinRange(mesh, idx, 3);
  Eigen::Matrix3d tmat = getVertexTransformationMatrix(mesh, idx);
  Eigen::Vector3d origin;
  origin << mesh.vertices[idx][0], mesh.vertices[idx][1], mesh.vertices[idx][2];
  Eigen::MatrixXd data, jacobian;
  data.resize(included_indices.size(), 3);
  data.resize(included_indices.size(), 2);
  for(int i = 0; i < included_indices.size(); i++){
    Eigen::Vector3d position;
    position << mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2];
    position = position - origin;
    //add shifted, rotated position to data vector
    data.row(i) = tmat*position;
  }
  Functor2Form minfunc(data);
  Eigen::LevenbergMarquardt<Functor2Form> lm_algo(minfunc);
  Eigen::VectorXd x;
  x << 1, 1, 1; 
  int info = lm_algo.minimize(x);

  eval << x(0), x(1), x(1), x(2);
  return eval;
}

//computes the curvature at each vertex by fitting the second fundamental form
//returns a vector of pairs containing the principal curvatures of a given vertex
std::vector<std::array<double,2> > computeMeshCurvature(const Mesh& mesh){
  std::vector<std::array<double,2> > curvatures;
  for(int i = 0; i < mesh.nvtx; i++){
    Eigen::Matrix2d eval = get2ndFormTensor(mesh, i);
    Eigen::EigenSolver<Eigen::Matrix2d> solver(eval);
    std::array<double,2> arr = {solver.eigenvalues().real()[0], solver.eigenvalues().real()[1]};
    curvatures.push_back(arr);
  }
  return curvatures; 
}

std::string printPLYWithCurvature(const Mesh& mesh){
  std::vector<std::array<double,2>> curvatures = computeMeshCurvature(mesh);
  double max = -1e10;
  for(int i = 0; i < curvatures.size(); i++){
    for(int j = 0; j < 2; j++){
      if(fabs(curvatures[i][j]) > max) max = curvatures[i][j];
    }
  }
  std::vector<double> gaussian_curvatures(curvatures.size());
  std::vector<double> average_curvatures(curvatures.size());
  for(int i = 0; i < curvatures.size(); i++){
    gaussian_curvatures[i] = curvatures[i][0] * curvatures[i][1] / (max*max);
    average_curvatures[i] = 0.5*(curvatures[i][0] + curvatures[i][1]) / max;
  }

  std::stringstream ss;
  ss << "ply\nformat ascii 1.0\n";
  ss << "element vertex " << mesh.nvtx << "\n";
  ss << "property float x\n";
  ss << "property float y\n";
  ss << "property float z\n";
  ss << "property uchar red\n";
  ss << "property uchar green\n";
  ss << "property uchar blue\n";
  ss << "element face " << mesh.ntri << "\n";
  ss << "property list uchar int vertex_index\n";
  ss << "end_header\n";
  for(int i = 0; i < mesh.nvtx; i++){
    ss << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << " ";
    double red_strength = 0.5*average_curvatures[i] + 0.5;
    double blue_strength = -0.5*average_curvatures[i] + 0.5;
    ss << round(255*red_strength) << " " << "0" << round(255*blue_strength) << "\n";
  }
  for(int i = 0; i < mesh.ntri; i++){
    ss << "3 " << mesh.triangles[i].indices[0] << " " << mesh.triangles[i].indices[1] << " " << mesh.triangles[i].indices[2] << "\n";
  }

  return ss.str();
}