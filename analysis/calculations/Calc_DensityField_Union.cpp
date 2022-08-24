#include "Calc_DensityField_Union.hpp"
#include "Eigen/Eigen"
void Calc_DensityField_Union::finalOutput(){
  for(int i = 0; i < df_vec_.size(); i++){
    FANCY_ASSERT(data_matrix_.size() == df_vec_[i]->size(), "Not all DensityField are of the appropriate size...");
  }
  data_matrix_.fill(0.0); 
  std::vector<DensityFieldDataPack> densityfields;
  for(auto df : df_vec_){
    densityfields.push_back(df->getAverageDatapack());
  }
  #pragma omp parallel for
  for(std::size_t i = 0; i < data_matrix_.size_1d(); i++){
    auto idx3d = data_matrix_.map1N(i);
    std::array<int,3> idx3d_int = {(int)idx3d[0], (int)idx3d[1], (int)idx3d[2]};
    for(auto& field : densityfields){
      data_matrix_.at(idx3d) += field.gridvals_[field._map31(idx3d_int)];
    }
  }
  //set the union of densityfields to the sum of each of the df values
  //now get the gradient at each point, will output its magnitude
  auto grad = gradient(data_matrix_);
  Matrix<double,3> gradient_magnitudes(data_matrix_.size());
  Matrix<double,3> hessian_magnitudes(data_matrix_.size()); //gaussian positive/negative definite
  Matrix<double,3> hessian_magnitudes3(data_matrix_.size()); //gaussian saddles ++-, +--
  Matrix<double,3> hessian_magnitudes2(data_matrix_.size()); //mean curvature

  for(std::size_t i = 0; i < data_matrix_.size_1d(); i++){
    auto gradient = grad.at_1d(i);
    double grad_step = 0.0;
    for(auto val : gradient){
      grad_step+= val*val;
    }
    grad_step = std::sqrt(grad_step);
    gradient_magnitudes.at_1d(i) = grad_step;
  }
  //now calculate the Hessian at each point, will need to calculate the Eigenvalues
  auto hessian = Hessian(data_matrix_);
  for(std::size_t i = 0; i < data_matrix_.size_1d(); i++){
    auto hessian_step = hessian.at_1d(i);
    Eigen::Matrix3d hessian_eigen;
    for(std::size_t i2 = 0; i2 < 3; i2++){
      for(std::size_t j = 0; j < 3; j++){
        std::array<std::size_t, 2> idx = {i2,j};
        hessian_eigen(i2,j) = hessian_step.at(idx);
      }
    }
    hessian_eigen = 0.5*(hessian_eigen.transpose()+hessian_eigen);
    Eigen::EigenSolver<Eigen::Matrix3d> es(hessian_eigen);
    auto ev = es.eigenvalues().real();
    double val = ev[0]*ev[1]*ev[2];
    double val3 = val;
    //zero all values that are neither positive definite or negative definite
    if(  !((ev[0] > 0.0 && ev[1] > 0.0 && ev[2] > 0.0) || (ev[0] < 0.0 && ev[1] < 0.0 && ev[2] < 0.0))  ){
      val = 0.0;
    }
    double val2 = ev.mean();
    hessian_magnitudes.at_1d(i) = val;//mean_eigenvalue;
    hessian_magnitudes2.at_1d(i) = val2;

    if( ((ev[0] > 0.0 && ev[1] > 0.0 && ev[2] > 0.0) || (ev[0] < 0.0 && ev[1] < 0.0 && ev[2] < 0.0))  ){
      val3 = 0.0;
    }
    hessian_magnitudes3.at_1d(i) = val3;

  }
  //now have 3 matrices, the sum of each of the fields, their gradient, and their hessian
  //need to output to separate files

  std::ofstream ofile(base_ + "_DensityField.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
  auto npoints = data_matrix_.size();
  ofile << "\%nx = [ " << npoints[0] << " " << npoints[1] << " " << npoints[2] << " ]\n";
  for(std::size_t k = 0; k < npoints[2]; k++){
  for(std::size_t j = 0; j < npoints[1]; j++){
  for(std::size_t i = 0; i < npoints[0]; i++){
        Vec3<std::size_t> pos = {i,j,k};
        ofile << data_matrix_.at(pos) << "   ";
      }
    }
  }
  ofile.close();
  ofile.open(base_ + "_GradientField.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
  ofile << "\%nx = [ " << npoints[0] << " " << npoints[1] << " " << npoints[2] << " ]\n";
  for(std::size_t k = 0; k < npoints[2]; k++){
  for(std::size_t j = 0; j < npoints[1]; j++){
  for(std::size_t i = 0; i < npoints[0]; i++){
        Vec3<std::size_t> pos = {i,j,k};
        ofile << gradient_magnitudes.at(pos) << "   ";
      }
    }
  }
  ofile.close();
  ofile.open(base_ + "_DefiniteGaussian.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
  ofile << "\%nx = [ " << npoints[0] << " " << npoints[1] << " " << npoints[2] << " ]\n";
  for(std::size_t k = 0; k < npoints[2]; k++){
  for(std::size_t j = 0; j < npoints[1]; j++){
  for(std::size_t i = 0; i < npoints[0]; i++){
        Vec3<std::size_t> pos = {i,j,k};
        ofile << hessian_magnitudes.at(pos) << "   ";
      }
    }
  }
  ofile.close();
  ofile.open(base_ + "_SaddleGaussian.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
  ofile << "\%nx = [ " << npoints[0] << " " << npoints[1] << " " << npoints[2] << " ]\n";
  for(std::size_t k = 0; k < npoints[2]; k++){
  for(std::size_t j = 0; j < npoints[1]; j++){
  for(std::size_t i = 0; i < npoints[0]; i++){
        Vec3<std::size_t> pos = {i,j,k};
        ofile << hessian_magnitudes3.at(pos) << "   ";
      }
    }
  }
  ofile.close();
  ofile.open(base_ + "_MeanCurvature.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
  ofile << "\%nx = [ " << npoints[0] << " " << npoints[1] << " " << npoints[2] << " ]\n";
  for(std::size_t k = 0; k < npoints[2]; k++){
  for(std::size_t j = 0; j < npoints[1]; j++){
  for(std::size_t i = 0; i < npoints[0]; i++){
        Vec3<std::size_t> pos = {i,j,k};
        ofile << hessian_magnitudes2.at(pos) << "   ";
      }
    }
  }
  return;
}