#include "Calc_SWIPES_CircleFit.hpp"
#include "Eigen/Eigen"
#include "unsupported/Eigen/NonLinearOptimization"
/*
//fits a circle to a 2d grid density
class Calc_CircleFit : public Calculation{
public:
  Calc_CircleFit(InputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
private:
  Vec<Vec<double> > angle_distribution_;
  Vec<double> angle_vec_;
  Vec<int> step_vec_;
  Vec<double> time_vec_;
  double mean_, var_, value_;
  //this one depeends on 3 atomgroups of equal size
  Vec3<AtomGroup*> atom_groups_;
};*/

struct FunctorCircleFit{
  Eigen::MatrixXd data;
  int values_;
  std::array<bool,3> fix_params_;
  FunctorCircleFit(const Eigen::MatrixXd& data_in, std::array<bool,3> fix_params){
    data = data_in;
    values_ = data_in.rows();
    fix_params_ = fix_params;
    return;
  }
  int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec){
    assert(b.size() == 3);
    assert(fvec.size() == values_);
    for(int i = 0; i < fvec.size(); i++){
      fvec(i) =  std::pow(data(i,0) - b(0), 2) + std::pow(data(i,1) - b(1),2)  - b(2);
    }
    return 0;
  }
  int df(const Eigen::VectorXd &b, Eigen::MatrixXd &fjac)
  {
    assert(b.size() == 3);
    assert(fjac.rows() == values_);
    for(int i = 0; i < data.rows(); i++){
      double dx = data(i,0) - b(0);
      double dy = data(i,1) - b(1);
      double j1 = -2.0*b(0)*dx;
      double j2 = -2.0*b(1)*dy;
      double j3 = -1.0;
      Eigen::Vector3d jac_row;
      jac_row << j1, j2, j3;
      fjac.row(i) = jac_row;
    }
    return 0;
  }
  int values(){return values_;}
};


Calc_SWIPES_CircleFit::Calc_SWIPES_CircleFit(InputPack& input) : Calculation{input} {
  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");

  std::string calc_name;
  input.params().readString("calculation", KeyType::Required, calc_name);
  calc_ = dynamic_cast<Calc_2D_Density*>(input.findCalculation(calc_name));
  FANCY_ASSERT(calc_ != 0, "Failed to find specified calculation.");

  Vec<double> xmin, xmax;
  input.params().readVector("xmin", KeyType::Required, xmin);
  input.params().readVector("xmax", KeyType::Required, xmax);
  FANCY_ASSERT(xmin.size() == 2, "Invalid xmin size in CircleFit object.");
  FANCY_ASSERT(xmax.size() == 2, "Invalid xmax size in CircleFit object.");
  for(int i = 0; i < 2; i++){
    xmin_[i] = xmin[i];
    xmax_[i] = xmax[i];
  }

  std::vector<double> params(3, 1.0);
  input.params().readVector("params", KeyType::Optional, params);
  FANCY_ASSERT(params.size() == 3, "Invalid params size in CircleFit object.");
  for(int i = 0; i < 3; i++){
    params_[i] = params[i];
  }

  std::vector<bool> fix_flags(3, 0);
  input.params().readVector("fix_params", KeyType::Optional, fix_flags);
  FANCY_ASSERT(fix_flags.size() == 3, "Invalid fix_params size in CircleFit object.");
  for(int i = 0; i < 3; i++){
    fix_params_[i] = fix_flags[i];
  } 
  input.params().readNumber("density", KeyType::Required, density_);
  input.params().readNumber("normal_direction", KeyType::Required, normal_direction_);
  FANCY_ASSERT(normal_direction_ == 1 || normal_direction_ == 0, "Invalid normal direction specified.");
  if(normal_direction_ == 1) parallel_direction_ = 0;
  else parallel_direction_ = 1;

  return;

}
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
void Calc_SWIPES_CircleFit::finalOutput(){
  auto g = calc_->getGridSubsection(xmin_, xmax_);
  //simple linear interpolation
  vertex_positions_.resize(g.grid_size[parallel_direction_], {0.0,0.0});
  if(normal_direction_ == 0){
    for(int j = 0; j < g.grid_size[1]; j++){
      double last_value = g.grid_data[0][j];
      double lzd = last_value - 0.5*density_;
      for(int i = 1; i < g.grid_size[0]; i++){
        double current_value = g.grid_data[i][j];
        double czd = current_value - 0.5*density_;
        //if the sign changes then we've crossed a threshold
        if(sgn(czd) != sgn(lzd)){
          //using linear interpolation, predict where the zero density is
          double x1 = (i-1)*g.grid_spacing[0];
          double x2 = i*g.grid_spacing[0];
          double y1 = lzd;
          double y2 = czd;
          double x = -y1*(x2-x1)/(y2-y1) + x1;
          double y = j*g.grid_spacing[1] + 0.5*g.grid_spacing[1];
          vertex_positions_[j] = {x + g.real_offset[0],y + g.real_offset[1]};
          break;
        }
        last_value = current_value;
        lzd = czd;
      }
    }
  }
  if(normal_direction_ == 1){
    for(int i = 0; i < g.grid_size[0]; i++){
      double last_value = g.grid_data[i][0];
      double lzd = last_value - 0.5*density_;
      bool crossed_threshold = 0;
      for(int j = 1; j < g.grid_size[1]; j++){
        double current_value = g.grid_data[i][j];
        double czd = current_value - 0.5*density_;
        //if the sign changes then we've crossed a threshold
        if(sgn(current_value-(0.5*density_)) != sgn(last_value-(0.5*density_))){
          //using linear interpolation, predict where the zero density is
          double x1 = (i-1)*g.grid_spacing[1];
          double x2 = i*g.grid_spacing[1];
          double y1 = lzd;
          double y2 = czd;
          double y = -y1*(x2-x1)/(y2-y1) + x1;
          double x = i*g.grid_spacing[0] + 0.5*g.grid_spacing[0];
          vertex_positions_[i] = {x + g.real_offset[0],y + g.real_offset[1]};
          break;
        }
        last_value = current_value;
        lzd = czd;
      }
    }
  }

  //vertex positions should correspond to an actual grid now, just have to set up fit
  Eigen::MatrixXd data(vertex_positions_.size(), 2);
  Eigen::VectorXd b(3);
  b << params_[0], params_[1], params_[2];
  for(int i = 0; i < vertex_positions_.size(); i++){
    data(i,0) = vertex_positions_[i][0];
    data(i,1) = vertex_positions_[i][1];
  }
  FunctorCircleFit minfunc(data, fix_params_);
  Eigen::LevenbergMarquardt<FunctorCircleFit> lm_algo(minfunc);
  int info = lm_algo.minimize(b);  
  for(int i = 0; i < 3; i++) params_[i] = b(i);
  //should be able to return circle data now




  //now I have a fitted circle and need to determine which side of the interface to cut off
  //get com of vertices
  double com = 0.0;
  for(int i = 0; i < vertex_positions_.size(); i++){
    com += vertex_positions_[i][normal_direction_];
  }
  com *= 1.0/vertex_positions_.size();

  //can be any number 1-4. 1-2 are axis1 aligned with 
  int direction = 0; //use the lower half of the sphere
  if(params_[normal_direction_] < com) direction = 1; //use the upper half of the sphere

  std::ofstream ofile(base_+"_circlefit.txt");
  //these parameters can be used to construct a new probe volume
  //above/below determines whether or not a union or intersection is appropriate
  ofile << "x0 = " << params_[0] << ",     y0 = " << params_[1] << ",     R = " << params_[2];
  ofile << ",     axis = " << g.integration_direction << ",     union = " << direction << std::endl;
  std::cout << "\n\n";
  for(int i = 0; i < vertex_positions_.size(); i++){
    ofile << vertex_positions_[i][0] << "     " << vertex_positions_[i][1] << std::endl;
  }
  ofile.close();

  return;
}