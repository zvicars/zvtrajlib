#include "Calc_1D_Density.hpp"
#include "Eigen/Eigen"
#include "unsupported/Eigen/NonLinearOptimization"

Calc_1D_Density::Calc_1D_Density(InputPack& input) : Calculation{input} {

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  input.params().readNumber("npoints", KeyType::Required, npoints_);
  input.params().readNumber("axis", KeyType::Required, dim_);
  fitSigmoidal = 0;
  input.params().readFlag("fit", KeyType::Optional, fitSigmoidal);
  frame_counter_ = 0;
  std::vector<int> xrange;
  input.params().readVector("x_range", KeyType::Optional, xrange);
  FANCY_ASSERT(xrange.size() == 2 && xrange[1] > xrange[0], "Invalid xrange parameter set in 1D Density calculation.");
  guess_.resize(3, 1.0);
  input.params().readVector("guess", KeyType::Optional, guess_);
  FANCY_ASSERT(xrange.size() == 2, "Invalid guess provided");
  idx_range_[0] = xrange[0];
  idx_range_[1] = xrange[1];
  frame_counter_ = 0;  
  grid_density_.resize(npoints_, 0.0);
  average_grid_density_.resize(npoints_, 0.0);
  params_.resize(3, 0.0);
  return;
}

struct FunctorSigmoidalFit{
  Eigen::MatrixXd data;
  int values_;
  FunctorSigmoidalFit(const Eigen::MatrixXd& data_in){
    data = data_in;
    values_ = data_in.rows();
    return;
  }
  int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec){
    assert(b.size() == 3);
    assert(fvec.size() == values_);
    for(int i = 0; i < fvec.size(); i++){
      fvec(i) = b(0)/(1+exp(b(1)*(data(i,0) - b(2)))) - data(i,1);
    }
    return 0;
  }
  int df(const Eigen::VectorXd &b, Eigen::MatrixXd &fjac)
  {
    assert(b.size() == 3);
    assert(fjac.rows() == values_);
    for(int i = 0; i < data.rows(); i++){
      Eigen::Vector3d jac_row;
      double dx = data(i,0) - b(2);
      jac_row << 1/(1+exp(b(1)*dx)),  -b(0)*(dx)*exp(b(1)*dx)/std::pow( 1 + exp(b(1)*dx), 2), b(0)*b(1)*exp(b(1)*dx)/std::pow( 1 + exp(b(1)*dx), 2);
      fjac.row(i) = jac_row;
    }
    return 0;
  }
  int values(){return values_;}
};

void Calc_1D_Density::calculate(){
  if(!doCalculate()) return;
  for(int i = 0; i < atom_group_->getIndices().size(); i++ ){
    int idx = atom_group_->getIndices()[i];
    auto position = box->atoms[idx].x;
    putInBin(position);
  }
  for(int i = 0; i < grid_density_.size(); i++){
    average_grid_density_[i] += grid_density_[i];
  }
  average_grid_spacing_ += grid_spacing_;
  if(fitSigmoidal){
  Eigen::MatrixXd data(idx_range_[1] - idx_range_[0] + 1, 2);
  int iterator = 0;
  for(int i = idx_range_[0]; i <= idx_range_[1]; i++){
    data(iterator, 1) = grid_density_[i];
    data(iterator, 0) = (iterator+0.5)*grid_spacing_;
    iterator++;
  }
  
  FunctorSigmoidalFit f1(data);
  Eigen::LevenbergMarquardt<FunctorSigmoidalFit> lm_algo(f1);
  Eigen::VectorXd b(3);
  b << guess_[0], guess_[1], guess_[2] - idx_range_[0]*grid_spacing_;
  int info = lm_algo.minimize(b);  
  b(2) += grid_spacing_*idx_range_[0];
  params_[0] = b(0);
  params_[1] = b(1);
  params_[2] = b(2);
  fits_.push_back(params_);
  tvec_.push_back(current_time_);
  frame_vec_.push_back(current_frame_);
  }
  frame_counter_++;
}
void Calc_1D_Density::update(){
  Calculation::update();
  for(int i = 0; i < grid_density_.size(); i++){
    grid_density_[i] = 0.0;
  }
  box_size_ = box->boxvec[dim_][dim_];
  grid_spacing_ = box_size_ / (double)(npoints_);
  return;
}
std::string Calc_1D_Density::printConsoleReport(){
  return "";
}
void Calc_1D_Density::finalOutput(){
  for(int i = 0; i < average_grid_density_.size(); i++){
    average_grid_density_[i] *= 1.0/(double)frame_counter_;
  }
  average_grid_spacing_*= 1.0/frame_counter_;
  if(fitSigmoidal){
  Eigen::MatrixXd data(idx_range_[1] - idx_range_[0] + 1, 2);
  int iterator = 0;
  for(int i = idx_range_[0]; i <= idx_range_[1]; i++){
    data(iterator, 1) = average_grid_density_[i];
    data(iterator, 0) = (iterator+0.5)*average_grid_spacing_;
    iterator++;
  }
  
  FunctorSigmoidalFit f1(data);
  Eigen::LevenbergMarquardt<FunctorSigmoidalFit> lm_algo(f1);
  Eigen::VectorXd b(3);
  b << guess_[0], guess_[1], guess_[2] - idx_range_[0]*average_grid_spacing_;
  int info = lm_algo.minimize(b);  
  b(2) += average_grid_spacing_*idx_range_[0];
  params_[0] = b(0);
  params_[1] = b(1);
  params_[2] = b(2);
  }
  std::ofstream ofile(base_ + "_avg_sigmoidal.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for 1D density calculation.");
  if(fitSigmoidal) ofile << "# " << params_[0] << "   " << params_[1] << "   " <<  params_[2] << std::endl;
  for(int i = 0; i < average_grid_density_.size(); i++){
    ofile << (i+0.5)*average_grid_spacing_ << "     " << average_grid_density_[i] << std::endl;
  }
  ofile.close();
  if(fitSigmoidal){
    ofile.open(base_ + "_ts_sigmoidal.txt");
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for 1D density calculation.");
    for(int i = 0; i < tvec_.size(); i++){
      ofile << tvec_[i] << "   " << frame_vec_[i] << "   ";
      for(int j = 0; j < 3; j++){
        ofile << fits_[i][j] << "   ";
      }
      ofile << "\n";
    }
    ofile.close();
  }
  return;
}