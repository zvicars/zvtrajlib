#include "PV_Cylinder.hpp"

Eigen::Matrix3d getRotationMatrix(Eigen::Vector3d cylinder_axis){
  auto normal_in = cylinder_axis;
  Eigen::Vector3d vertex, normal;
  Eigen::Vector3d unit_z;
  unit_z << 0, 0, 1;
  for(int i= 0; i < 3; i++){
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

PV_Cylinder::PV_Cylinder(InputPack& input):ProbeVolume{input}
{
  Vec<double> position, axis;
  input.params().readVector("position", KeyType::Required, position);
  input.params().readVector("axis", KeyType::Required, axis);
  input.params().readNumber("radius", KeyType::Required, radius_);
  input.params().readNumber("height", KeyType::Required, height_);

  FANCY_ASSERT(position.size() == 3, "Invalid position vector size in ProbeVolume " + name_);
  FANCY_ASSERT(axis.size() == 3, "Invalid axis vector size in ProbeVolume " + name_);

  for(int i = 0; i < 3; i++){
    axis_(i) = axis[i];
    position_(i) = position[i];
  }
  tmat_ = getRotationMatrix(axis_);
  return;
}

double PV_Cylinder::compute(Vec3<double> position){
  Eigen::Vector3d x, x_trans;
  x << position[0], position[1], position[2];
  x_trans = tmat_*(x - position_);
  
  double dist =sqrt(x_trans(0)*x_trans(0) + x_trans(1)*x_trans(1));
  if(dist <= radius_ && x_trans(2) >= 0 && x_trans(2) <= height_) return 1.0;
  return 0.0;
}