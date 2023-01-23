#include "Calc_1D_Density_IP.hpp"
#include "../helper/functors.hpp"
#include "../../tools/pbcfunctions.hpp"
inline double heaviside(double x){
	if(x <= 0) return 0;
	return 1.0; 
}

inline double h_x(double x, double xmin, double xmax, double sigma, double xc){
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


Calc_1D_Density_IP::Calc_1D_Density_IP(InputPack& input) : Calculation{input} {

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  
  input.params().readNumber("npoints", KeyType::Required, npoints_);
  input.params().readNumber("axis", KeyType::Required, dim_);
  fitSigmoidal = 0;
  input.params().readFlag("fit", KeyType::Optional, fitSigmoidal);

  std::vector<int> xrange;
  input.params().readVector("x_range", KeyType::Optional, xrange);
  FANCY_ASSERT(xrange.size() == 2 && xrange[1] > xrange[0], "Invalid xrange parameter set in 1D Density calculation.");

  input.params().readVector("guess", KeyType::Optional, guess_);
  FANCY_ASSERT(guess_.size() == 5, "Invalid guess provided");
  fix_.fill(0);
  std::vector<bool> fix;
  input.params().readVector("fix", KeyType::Optional, fix);
  FANCY_ASSERT(fix_.size() == 5, "Invalid guess provided");  
  for(int i = 0; i < fix.size(); i++){
    fix_[i] = fix[i];
  }
  idx_range_[0] = xrange[0];
  idx_range_[1] = xrange[1];
  
  //com_corr determines center-of-mass correction behavior
  //after computing a 1d density field, it will translate it such that the maximum density
  //is either in the middle of the box 1, the left-edge of the box 2 or it stays uncorrected 0
  com_corr_ = 0;
  input.params().readNumber("comcorrect", KeyType::Optional, com_corr_);

  frame_counter_ = 0;
  grid_density_.resize(npoints_, 0.0);
  average_grid_density_.resize(npoints_, 0.0);

  for(int i = 0; i < 5; i++){
    params_[i] = guess_[i];
  }

  coarseGrain = 0;
  input.params().readFlag("smear", KeyType::Optional, coarseGrain);
  if(coarseGrain){
    input.params().readNumber("sigma", KeyType::Required, sigma_);
  }
  std::string pv_name;
  input.params().readString("probe_volume", KeyType::Required, pv_name);
  auto pv_pointer = input.findProbeVolume(pv_name);
  FANCY_ASSERT(pv_pointer != 0, "Failed to find specified probe volume.");
  pv_ = pv_pointer;

  return;
}


void Calc_1D_Density_IP::add_gaussian(double x_in)
{
    double x = x_in;
    int lxmin = floor((x-2*sigma_)/grid_spacing_);
    int lxmax = ceil((x+2*sigma_)/grid_spacing_); 
    #pragma omp parallel for
    for(int ix = lxmin; ix <= lxmax; ix++)
    {
      int idx = ix;
      if(idx >= npoints_) idx -= npoints_;
      else if(idx < 0) idx += npoints_;
      double xmin, xmax;
      xmin = ix * grid_spacing_;
      xmax = xmin + grid_spacing_;
      grid_density_[idx] += h_x(x, xmin, xmax, sigma_, 2.0*sigma_);
    }
    return;
}

void Calc_1D_Density_IP::calculate(){
  if(!doCalculate()){
    return;
  }
  auto indices = atom_group_->getIndices();
  int counter = 0;
  double com_dx_ = 0.0;
  if(com_corr_ != 0){
    //first, figure out if split near edge of box by comparing variances for 
    //normal and half-box shifted configurations
    double var0=0.0, var1=0.0, com0=0.0, com1=0.0;
    for(auto idx : indices){
      if(pv_->compute(box->atoms[idx].x)==0.0) continue;
      com0 += wrapNumber(box->atoms[idx].x[dim_], box_size_);
      com1 += wrapNumber(box->atoms[idx].x[dim_] + 0.5*box_size_, box_size_);
    }
    com0 /= (double)indices.size();
    com1 /= (double)indices.size();
    //compute mean-centered variances
    for(auto idx : indices){
      if(pv_->compute(box->atoms[idx].x)==0.0) continue;
      var0 += pow(wrapNumber(box->atoms[idx].x[dim_], box_size_) - com0, 2);
      var1 += pow(wrapNumber(box->atoms[idx].x[dim_] + 0.5*box_size_, box_size_) - com1, 2);
    }
    var0 /= (double)indices.size();
    var1 /= (double)indices.size();
    //if half-length shifting gives lower variance, then we'll want to bake in an initial half-length shift
    //suggests that density is straddling box edge
    if(var1 < var0) com_dx_ = 0.5*box_size_ - com1;
    else com_dx_ = -com0;
    //if comm_corr_ == 1, we want the high density to be in the center of the box instead of 0
    if(com_corr_ == 1) com_dx_ += 0.5*box_size_;
    com_dx_ = com_dx_;
  }
  for(auto idx : indices){
    counter++;
    if(pv_->compute(box->atoms[idx].x)==0.0) continue;
    bool out_of_range_flag = 1;
    if(idx >= 0 && idx < box->atoms.size()) out_of_range_flag = 0;
    if(out_of_range_flag == 1){
      std::ofstream ofile("dump_atomgroup_info.txt");
      ofile << atom_group_ << std::endl;
      std::string ag_dump = atom_group_->getDumpString();
      for(auto idx2 : indices) ofile << idx2 << "   " << "\n";
      ofile << ag_dump;
      ofile.close();
    }
    FANCY_ASSERT(!out_of_range_flag, "AtomGroup provided index "
     + std::to_string(idx)
     + " which is not compatible with the number of atoms in the box. This was done on iteration " + std::to_string(counter) + " at time "
     + std::to_string(box->time));

    auto position = box->atoms[idx].x;
    position[dim_] += com_dx_;
    position[dim_] = wrapNumber(position[dim_], box_size_);
    if(coarseGrain) add_gaussian(position[dim_]);
    else putInBin(position);
  }

  for(int i = 0; i < grid_density_.size(); i++){
    average_grid_density_[i] += grid_density_[i];
  }
  average_grid_spacing_ += grid_spacing_;
  performFitStep();
  frame_counter_++;
  return;
}
void Calc_1D_Density_IP::update(){
  if(hasUpdated()) return;
  Calculation::update();
  for(auto& entry : grid_density_){
    entry = 0.0;
  }
  box_size_ = box->boxvec[dim_][dim_];
  grid_spacing_ = box_size_ / (double)(npoints_);
  return;
}
std::string Calc_1D_Density_IP::printConsoleReport(){
  return "";
}
void Calc_1D_Density_IP::finalOutput(){
  for(auto& gridval : average_grid_density_){
    gridval *= 1.0/(double)frame_counter_;
  }
  average_grid_spacing_*= 1.0/frame_counter_;
  performFitAvg();
  std::ofstream ofile(base_ + "_avg_sigmoidal.txt");
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for 1D density calculation.");
  if(fitSigmoidal) ofile << "# " << params_[0] << "   " << params_[1] << "   " <<  params_[2] << "   " 
                         <<  params_[3] << "   " << params_[4] << std::endl;
  for(int i = 0; i < average_grid_density_.size(); i++){
    ofile << (i+0.5)*average_grid_spacing_ << "     " << average_grid_density_[i] << std::endl;
  }
  ofile.close();

  if(fitSigmoidal){
    ofile.open(base_ + "_ts_sigmoidal.txt");
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for 1D density calculation.");
    for(int i = 0; i < fits_.size(); i++){
      ofile << tvec_[i] << "   " << frame_vec_[i] << "   ";
      for(int j = 0; j < 5; j++){
        ofile << fits_[i][j] << "   ";
      }
      ofile << "\n";
    }
    ofile.close();
  }
  return;
}

void Calc_1D_Density_IP::performFitStep(){
  if(!fitSigmoidal) return;
  Eigen::MatrixXd data(idx_range_[1] - idx_range_[0] + 1, 2);
  int iterator = 0;
  for(int i = idx_range_[0]; i <= idx_range_[1]; i++){
    data(iterator, 1) = grid_density_[i];
    data(iterator, 0) = (iterator+idx_range_[0]+0.5)*grid_spacing_;
    iterator++;
  }
  logisticStepFunctor f1(data, fix_);
  Eigen::LevenbergMarquardt<logisticStepFunctor> lm_algo(f1);
  Eigen::VectorXd b(5);
  b << guess_[0], guess_[1], guess_[2], guess_[3], guess_[4];
  int info = lm_algo.minimize(b);  
  params_[0] = b(0);
  params_[1] = b(1);
  params_[2] = b(2);
  params_[3] = b(3);
  params_[4] = b(4);
  fits_.push_back(params_);
  tvec_.push_back(current_time_);
  frame_vec_.push_back(current_frame_);
}

void Calc_1D_Density_IP::performFitAvg(){
  if(!fitSigmoidal) return;
  Eigen::MatrixXd data(idx_range_[1] - idx_range_[0] + 1, 2);
  int iterator = 0;
  for(int i = idx_range_[0]; i <= idx_range_[1]; i++){
    data(iterator, 1) = average_grid_density_[i];
    data(iterator, 0) = (iterator+idx_range_[0]+0.5)*average_grid_spacing_;
    iterator++;
  }
  
  logisticStepFunctor f1(data, fix_);
  Eigen::LevenbergMarquardt<logisticStepFunctor> lm_algo(f1);
  Eigen::VectorXd b(5);
  b << guess_[0], guess_[1], guess_[2], guess_[3], guess_[4];
  int info = lm_algo.minimize(b);  
  params_[0] = b(0);
  params_[1] = b(1);
  params_[2] = b(2);
  params_[3] = b(3);
  params_[4] = b(4);
  return;
}