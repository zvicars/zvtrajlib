#include "Calc_1D_DensityFull.hpp"
#include "../helper/functors.hpp"
#include "../../tools/stlmath.hpp"
inline double Calc_1D_DensityFull::heaviside(double x){
	if(x <= 0) return 0;
	return 1.0; 
}

double Calc_1D_DensityFull::h_x(double x, double xmin, double xmax, double sigma, double xc){
    double eval = 0.0;
    double k, k1, k2, invk;
    double sigma2 = sigma*sigma;
    k = sqrt(2*M_PI)*sigma*erf(xc/(sqrt(2)*sigma)) - 2*xc*exp(-(xc*xc)/(2*sigma2));
    invk = 1/k;
    k1 = invk*sqrt(0.5*M_PI*sigma2);
    k2 = invk*exp(-0.5*(xc*xc)/sigma2);
    eval = (k1 * erf((xmax-x)/(sqrt(2)*sigma)) - k2*(xmax-x) - 0.5)*heaviside(xc - fabs(xmax-x))
    + (k1 * erf((x-xmin)/(sqrt(2)*sigma)) - k2*(x-xmin) - 0.5)*heaviside(xc - fabs(x-xmin))
    + heaviside(xc + 0.5*(xmax-xmin) - fabs(x - 0.5*(xmin+xmax)));
    return eval;
}


Calc_1D_DensityFull::Calc_1D_DensityFull(InputPack& input) : Calculation{input} {

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
  hasRange_ = input.params().readVector("x_range", KeyType::Optional, xrange);
  if(hasRange_){
    FANCY_ASSERT(xrange.size() == 2 && xrange[1] > xrange[0], "Invalid xrange parameter set in 1D Density calculation.");
    FANCY_ASSERT(xrange.size() == 2, "Invalid guess provided");
    idx_range_[0] = xrange[0];
    idx_range_[1] = xrange[1];
  }
  else{
    idx_range_[0] = 0;
    idx_range_[1] = npoints_ - 1;
  }
  guess_.resize(4, 1.0);
  input.params().readVector("guess", KeyType::Optional, guess_);
  FANCY_ASSERT(guess_.size() == 4, "guess vector must be of size 4");
  frame_counter_ = 0;  
  grid_density_.resize(npoints_, 0.0);
  average_grid_density_.resize(npoints_, 0.0);
  params_.resize(4, 0.0);

  coarseGrain = 0;
  input.params().readFlag("smear", KeyType::Optional, coarseGrain);
  if(coarseGrain){
    input.params().readNumber("sigma", KeyType::Optional, sigma_);
  }

  std::string pv_name;
  hasPV_ = input.params().readString("probevolume", KeyType::Optional, pv_name);
  if(hasPV_){
    pv_ = input.findProbeVolume(pv_name);
    FANCY_ASSERT(pv_ != 0, "failed to find probe volume");
  }

  return;
}

void Calc_1D_DensityFull::add_gaussian(double x_in)
{
    double x = x_in;
    int lxmin = floor((x-2*sigma_)/grid_spacing_);
    int lxmax = ceil((x+2*sigma_)/grid_spacing_); 
    #pragma omp parallel for
    for(int ix = lxmin; ix <= lxmax; ix++)
    {
      int idx = ix;
      if(idx >= grid_density_.size()) continue;
      else if(idx < 0) continue;
      double xmin, xmax;
      xmin = ix * grid_spacing_;
      xmax = xmin + grid_spacing_;
      grid_density_[idx] += h_x(x, xmin, xmax, sigma_, 2.0*sigma_);
    }
    return;
}

void Calc_1D_DensityFull::calculate(){
  if(!doCalculate()) return;
  auto& indices = atom_group_->getIndices();
  for(auto idx : indices){
    auto position = box->atoms[idx].x;
      if(hasPV_){
        if(pv_->compute(position) ==0){
          continue;
        }
      }
    if(coarseGrain) add_gaussian(position[dim_]);
    else putInBin(position);
  }
  average_grid_density_ = average_grid_density_ + grid_density_;
  average_grid_spacing_ += grid_spacing_;

  if(fitSigmoidal){
  Eigen::MatrixXd data(idx_range_[1] - idx_range_[0] + 1, 2);
  int iterator = 0;
  for(int i = idx_range_[0]; i <= idx_range_[1]; i++){
    data(iterator, 1) = grid_density_[i];
    data(iterator, 0) = (iterator + idx_range_[0] + 0.5)*grid_spacing_;
    iterator++;
  }
  logisticFunctorFull f1(data);
  Eigen::LevenbergMarquardt<logisticFunctorFull> lm_algo(f1);
  Eigen::VectorXd b(4);
  b << guess_[0], guess_[1], guess_[2], guess_[3];
  int info = lm_algo.minimize(b);  
  params_[0] = b(0);
  params_[1] = b(1);
  params_[2] = b(2);
  params_[3] = b(3);
  fits_.push_back(params_);
  tvec_.push_back(current_time_);
  frame_vec_.push_back(current_frame_);
  }

  frame_counter_++;
}
void Calc_1D_DensityFull::update(){
  if(hasUpdated()) return;
  Calculation::update();
  for(auto& val : grid_density_){
    val = 0.0;
  }
  box_size_ = box->boxvec[dim_][dim_];
  grid_spacing_ = box_size_ / (double)(npoints_);
  return;
}
std::string Calc_1D_DensityFull::printConsoleReport(){
  return "";
}
void Calc_1D_DensityFull::finalOutput(){
  for(int i = 0; i < average_grid_density_.size(); i++){
    average_grid_density_[i] *= 1.0/(double)frame_counter_;
  }
  average_grid_spacing_*= 1.0/frame_counter_;
  if(fitSigmoidal){
  Eigen::MatrixXd data(idx_range_[1] - idx_range_[0] + 1, 2);
  int iterator = 0;
  for(int i = idx_range_[0]; i <= idx_range_[1]; i++){
    data(iterator, 1) = average_grid_density_[i];
    data(iterator, 0) = (iterator+idx_range_[0]+0.5)*average_grid_spacing_;
    iterator++;
  }
  
  logisticFunctorFull f1(data);
  Eigen::LevenbergMarquardt<logisticFunctorFull> lm_algo(f1);
  Eigen::VectorXd b(4);
  b << guess_[0], guess_[1], guess_[2], guess_[3];
  int info = lm_algo.minimize(b);  
  params_[0] = b(0);
  params_[1] = b(1);
  params_[2] = b(2);
  params_[3] = b(3);
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
      for(int j = 0; j < 4; j++){
        ofile << fits_[i][j] << "   ";
      }
      ofile << "\n";
    }
    ofile.close();
  }
  return;
}