#include "MarchingCubesInterface.hpp"
#include "Eigen/Eigen"
#include "unsupported/Eigen/NonLinearOptimization"
#include <fstream>
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
  if(fabs(fabs(c) - 1) < 1e-5){
    rotation_matrix = Eigen::Matrix3d::Identity();
    if(c < 0){
      rotation_matrix(1,1) = -1;
      rotation_matrix(2,2) = -1;
    }
    return rotation_matrix;
  }
  kmat << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
  rotation_matrix = Eigen::Matrix3d::Identity() + kmat + kmat*kmat*((1-c)/(s*s));
  return rotation_matrix;
}

std::vector<std::vector<int> > buildNeighborList(const Mesh& mesh){
  std::vector<std::vector<int> > eval(mesh.nvtx);
  for(int idx = 0; idx < mesh.nvtx; idx++){
    std::vector<int> indices;
    for(int i = 0; i < mesh.triangles.size(); i++){
      for(int j = 0; j < 3; j++){
        if(mesh.triangles[i].indices[j] == idx){
          if(j == 0){
            indices.push_back(mesh.triangles[i].indices[1]);
            indices.push_back(mesh.triangles[i].indices[2]);
          }
          else if(j == 1){
            indices.push_back(mesh.triangles[i].indices[0]);
            indices.push_back(mesh.triangles[i].indices[2]);
          }
          else{
            indices.push_back(mesh.triangles[i].indices[0]);
            indices.push_back(mesh.triangles[i].indices[1]); 
          }
          break;
        }
      }
    }
    for(int i = 0; i < indices.size(); i++){
      for(int j = i+1; j < indices.size(); j++){
        if(indices[i] == indices[j]){
          indices.erase(indices.begin() + j);
          j--;
        }
      }
    }
    eval[idx] = indices;
  }
  return eval;
}

inline std::vector<int> getNeighborsWithinRange(const Mesh& mesh, const std::vector<std::vector<int> >& neighbor_list, int idx, int neighbors){
  std::vector<int> indices(1, idx);
  std::vector<int> last_indices(1, idx);
  std::vector<int> new_indices;
  for(int i = 0; i < neighbors; i++){
    for(int j = 0; j < last_indices.size(); j++){
      std::vector<int> new_index_step = neighbor_list[last_indices[j]];
      new_indices.insert(new_indices.end(), new_index_step.begin(), new_index_step.end());
    }
    for(int j = new_indices.size()-1; j >= 0; j--){
      for(int k = 0; k < indices.size(); k++){
        if(new_indices[j] == indices[k]){
          new_indices.erase(new_indices.begin() + j);
          break;
        }
      }
    }
    indices.insert(indices.end(), new_indices.begin(), new_indices.end());
    last_indices = new_indices;
    new_indices.clear();
  }
  return indices;
}

struct Functor2Form{
  Eigen::MatrixXd data;
  int values_;
  double invs_;
  Functor2Form(const Eigen::MatrixXd& data_in, double sens){
    data = data_in;
    values_ = data_in.rows();
    invs_ = 1.0/sens;
    return;
  }
  int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec){
    assert(b.size() == 3);
    assert(fvec.size() == values_);
    for(int i = 0; i < fvec.size(); i++){
      fvec(i) = invs_*(0.5*b(0)*data(i,0)*data(i,0) + b(1)*data(i,0)*data(i,1) + 0.5*b(2)*data(i,1)*data(i,1) - data(i,2));
    }
    return 0;
  }
  int df(const Eigen::VectorXd &b, Eigen::MatrixXd &fjac)
  {
    assert(b.size() == 3);
    assert(fjac.rows() == values_);
    for(int i = 0; i < data.rows(); i++){
      Eigen::Vector3d jac_row;
      jac_row << invs_*0.5*data(i,0)*data(i,0), invs_*data(i,0)*data(i,1), invs_*0.5*data(i,1)*data(i,1);
      fjac.row(i) = jac_row;
    }
    return 0;
  }
  int values(){return values_;}
};

Eigen::Matrix2d get2ndFormTensor(const Mesh& mesh, int idx, int neighbors, const std::vector<std::vector<int> >& nlist){
  Eigen::Matrix2d eval;
  std::vector<int> included_indices = getNeighborsWithinRange(mesh, nlist, idx, neighbors);
  Eigen::Matrix3d tmat = getVertexTransformationMatrix(mesh, idx);
  Eigen::Vector3d origin;
  origin << mesh.vertices[idx][0], mesh.vertices[idx][1], mesh.vertices[idx][2];
  Eigen::MatrixXd data;
  data.resize(included_indices.size(), 3);
  for(int i = 0; i < included_indices.size(); i++){
    Eigen::Vector3d position;
    int ii = included_indices[i];
    position << mesh.vertices[ii][0], mesh.vertices[ii][1], mesh.vertices[ii][2];
    position = position - origin;
    //add shifted, rotated position to data vector
    data.row(i) = tmat*position;
  }
  Functor2Form minfunc(data, 1.0);
  Eigen::LevenbergMarquardt<Functor2Form> lm_algo(minfunc);
  Eigen::VectorXd x;
  x.resize(3);
  x(0) = 0; x(1) = 0; x(2) = 0;
  int info = lm_algo.minimize(x);

  eval << x(0), x(1), x(1), x(2);
  /*
  if(idx%10 == 0){
    std::ofstream ofile("dump" + std::to_string(idx) + ".txt");
    ofile << "origin: " << origin;
    ofile << "vertex list:\n";
    for(int i = 0; i < included_indices.size(); i++){
      ofile << included_indices[i] << " ";
    }
    ofile << "\n";
    ofile << "fit results:\n" << x << "\n";
    ofile << "niter = " << lm_algo.iter << std::endl;

    for(int i = 0; i < included_indices.size(); i++){
      ofile << data.row(i) << "\n";
    }

    ofile.close();
  }
  */
  return eval;
}

//computes the curvature at each vertex by fitting the second fundamental form
//returns a vector of pairs containing the principal curvatures of a given vertex
std::vector<std::array<double,2> > computeMeshCurvature(const Mesh& mesh, int neighbors){
  std::vector<std::array<double,2> > curvatures;
  auto nlist = buildNeighborList(mesh);
  for(int i = 0; i < mesh.nvtx; i++){
    Eigen::Matrix2d eval = get2ndFormTensor(mesh, i, neighbors, nlist);
    Eigen::EigenSolver<Eigen::Matrix2d> solver(eval);
    std::array<double,2> arr = {solver.eigenvalues().real()[0], solver.eigenvalues().real()[1]};
    curvatures.push_back(arr);
  }
  return curvatures; 
}

std::string printPLYWithCurvature(const Mesh& mesh, const std::vector<double>& cv){
  double max = -1e10;
  for(int i = 0; i < cv.size(); i++){
      if(fabs(cv[i]) > max) max = fabs(cv[i]);
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
    double red_channel, green_channel;
    red_channel = 0.5*(cv[i]/max) + 0.5;
    green_channel = -0.5*(cv[i]/max) + 0.5;
    ss << round(255*red_channel) << " " << round(255*green_channel) << " 0" << "\n";
  }
  for(int i = 0; i < mesh.ntri; i++){
    ss << "3 " << mesh.triangles[i].indices[0] << " " << mesh.triangles[i].indices[1] << " " << mesh.triangles[i].indices[2] << "\n";
  }

  return ss.str();
}

std::string printPLYWithRGB(const Mesh& mesh, const std::vector<RGB>& cv){
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
    std::array<unsigned char,3> channels;
    channels[0] = cv[i].r;
    channels[1] = cv[i].g;
    channels[2] = cv[i].b;
    for(int j = 0; j < 3; j++){
      ss << (int)channels[j] << " ";
    }
    ss << "\n";
  }
  for(int i = 0; i < mesh.ntri; i++){
    ss << "3 " << mesh.triangles[i].indices[0] << " " << mesh.triangles[i].indices[1] << " " << mesh.triangles[i].indices[2] << "\n";
  }

  return ss.str();
}