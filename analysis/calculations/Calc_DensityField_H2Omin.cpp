#include "Calc_DensityField_H2Omin.hpp"
#include "../helper/functors.hpp"
#include "../../extern/LBFGS/LBFGS.h"
void Calc_DensityFieldH2OMin::calculate(){
  if(!doCalculate()) return;
  for(auto& value : gridvals_){
    value = 0.0;
  }
  for(int i = 0; i < 3; i++){
    span_[i] = ceil(cutoff_ / gridspacing_[i]);
  }
  auto atomIndices = atom_group_->getIndices();
  CellGrid c1(cutoff_, box_size_, 1);
  for(auto index : atomIndices){
    c1.addIndexToGrid(index, box->atoms[index].x);
  }
  #pragma omp parallel for
  for(int i = 0; i < gridvals_.size(); i++){
    auto index = _map13(i);
    Vec3<double> real_pos;
    for(int j = 0; j < 3; j++){
      real_pos[j] = (index[j]+0.5)*gridspacing_[j] + minx_[j];
    }
    auto& indices_step = c1.getNearbyIndices(real_pos);
    if(indices_step.size() == 0){
      gridvals_[i] = 0.0;
      continue;
    }
    Vec3<double> normal_step;
    gridvals_[i] = minimize_ff(real_pos, indices_step, normal_step);
    normals_[i] = normal_step;
  }
  nframes_++;
  avggridvals_ = avggridvals_ + gridvals_;
  avgNormals_ = avgNormals_ + normals_; 
  avggridspacing_ = avggridspacing_ + gridspacing_;
  return;
}


struct RotationMatrixPotentialFunctor
{
    Eigen::VectorXd x, y, z, q, e, s;
    Eigen::Vector3d pos;
    Eigen::Matrix3d rotationMatrix;
    double cutoff, beta, h_charge = 0.5897, m_charge = -2.0*0.5897;
    //initialized with parameter matrix
    //optimizes euler angles
    RotationMatrixPotentialFunctor(Eigen::MatrixXd data, Eigen::Vector3d pos_, double _cutoff, double _beta){
      //parameters are x, y, z, charge, epsilon, sigma, cutoff
      x = data.col(0);
      y = data.col(1);
      z = data.col(2);
      q = data.col(3);
      e = data.col(4);
      s = data.col(5);
      cutoff = _cutoff;
      beta = _beta;
      pos = pos_;
      return;
    };
    std::array<double,3> getNormal(Eigen::VectorXd &b){
      Eigen::Vector3d axis;
      setRotationMatrix(b);
      axis << -1,0,0;
      axis = rotationMatrix*axis;
      std::array<double,3> retVal;
      for(int i = 0; i < 3; i++){
        retVal[i] = axis[i];
      }
      return retVal;
    }
    //minimizing potential
    //first arg is x vector, second arg is the gradient, return value is potential
    double operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec)
    {
        assert(b.size()==3);
        setRotationMatrix(b);
        double retVal = f(b);
        df(b, fvec);
        return retVal;
    }
    void setRotationMatrix(const Eigen::VectorXd &b){
      Eigen::Vector3d radangles;
      radangles << b[2]*M_PI/180.0, b[1] * M_PI / 180.0, b[0] * M_PI / 180.0;
      Eigen::AngleAxisd rollAngle(radangles[0], Eigen::Vector3d::UnitZ());
      Eigen::AngleAxisd yawAngle(radangles[1], Eigen::Vector3d::UnitY());
      Eigen::AngleAxisd pitchAngle(radangles[2], Eigen::Vector3d::UnitX());
      Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;
      rotationMatrix = q.matrix(); 
      return;
    }
    double f(const Eigen::VectorXd &b){
      double eval_temp = 0.0; 
      Eigen::Vector3d x_step;
      Eigen::Vector3d pos_step = pos;
      for(int i = 0; i < x.size(); i++){
        x_step << x[i], y[i], z[i];
        pos_step = pos;
        pos_step(0) += 0.015; 
        rotateAboutPoint(pos_step, pos);
        auto r = (x_step - pos_step).norm();
        if(r < 0.2) r = 0.2;
        eval_temp += m_charge*(q(i)/r)*erfc(beta*r)*(138.935458*0.5);
        
        pos_step = pos;
        pos_step(0) += 0.6120792*0.09572; pos_step(1) += 0.790796*0.09572;
        rotateAboutPoint(pos_step, pos);
        r = (x_step - pos_step).norm();
        if(r < 0.2) r = 0.2;
        eval_temp += h_charge*(q(i)/r)*erfc(beta*r)*(138.935458*0.5);
        
        pos_step = pos;
        pos_step(0) += 0.6120792*0.09572; pos_step(1) -= 0.790796*0.09572;
        rotateAboutPoint(pos_step, pos);
        r = (x_step - pos_step).norm();
        if(r < 0.2) r = 0.2;
        eval_temp += h_charge*(q(i)/r)*erfc(beta*r)*(138.935458*0.5);

        pos_step = pos;
        r = (x_step - pos_step).norm();
        if(r < 0.2) r = 0.2;
        double r6 = std::pow(s(i)/r,6.0);
        double r12 = r6*r6;
        eval_temp += -4.0*e(i)*(r6 - r12);
      }
      return eval_temp;
    }
    void df(const Eigen::VectorXd &b, Eigen::VectorXd &fgrad){
      double dx = 0.01, inv_dx = 100.0;
      auto b_temp = b;
      double y0 = f(b_temp);
      for(int i = 0; i < 3; i++){
        b_temp(i) += dx;
        setRotationMatrix(b_temp);
        fgrad(i) = (f(b_temp) - y0) * inv_dx;
        b_temp(i) -= dx;
      }
      return;
    }
    void rotateAboutPoint(Eigen::Vector3d& x_step, const Eigen::Vector3d& x){
      x_step = rotationMatrix*(x_step-x) + x;
      return;
    }
};



double Calc_DensityFieldH2OMin::minimize_ff(Vec3<double> site_position, const std::vector<int>& indices, Vec3<double>& normal){
  //energy minimization procedure to get the optimal orientation
  Eigen::Vector3d pos;
  pos << site_position[0], site_position[1], site_position[2];
  Eigen::MatrixXd data;
  data.resize(indices.size(), 6);
  for(int i = 0; i < (int)indices.size(); i++){
    auto x_step = box->atoms[indices[i]].x;
    auto atom_step = box->atoms[indices[i]];
    data(i,0) = x_step[0];
    data(i,1) = x_step[1];
    data(i,2) = x_step[2];  
    data(i,3) = charge_map_.find(atom_step.name)->second;
    data(i,4) = eps_map_.find(atom_step.name)->second;
    data(i,5) = sigma_map_.find(atom_step.name)->second;
  }
  RotationMatrixPotentialFunctor fnct1(data, pos, cutoff_, beta_);
  LBFGSpp::LBFGSParam<double> param_;
  LBFGSpp::LBFGSSolver<double>* solver_;
  param_.epsilon = 1e-5;
  param_.max_iterations = 20;
  solver_ = new LBFGSpp::LBFGSSolver<double>(param_);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(3);
  double fx;
  int niter = 0;
  try{niter = solver_->minimize(fnct1, x, fx);}
  catch( ... ){};
  delete solver_;
  normal = fnct1.getNormal(x);
  return fx;
}
void Calc_DensityFieldH2OMin::finalOutput()
{
  Calc_DensityField::finalOutput();
  for(int i = 0; i < avgNormals_.size(); i++){
    avgNormals_[i] = (1.0/(double)nframes_) * avgNormals_[i];
  }
  std::ofstream ofile(base_ + "_normals.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for " + name_);
  for(int i = 0; i < npoints_[0]; i++){
  for(int j = 0; j < npoints_[1]; j++){
  for(int k = 0; k < npoints_[2]; k++){
        Vec3<int> pos = {i,j,k};
        int idx = _map31(pos);
        ofile << i << "   " << j << "   " << k << "   " << avgNormals_[idx][0] << "   " << avgNormals_[idx][1] 
        << "   " << avgNormals_[idx][2] << "\n";
      }
    }
  }
  ofile.close();  
  return;
}