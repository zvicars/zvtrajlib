#include "PV_INDUSRect.hpp"
#include "../../tools/smearfunctions.hpp"
PV_INDUSRect::PV_INDUSRect(InputPack& input):ProbeVolume{input}
{
  Vec<double> x_range, y_range, z_range;
  input.params().readVector("x_range", KeyType::Required, x_range);
  input.params().readVector("y_range", KeyType::Required, y_range);
  input.params().readVector("z_range", KeyType::Required, z_range);
  input.params().readNumber("sigma", KeyType::Required, sigma_);
  cutoff_ = 2.0*sigma_;
  input.params().readNumber("cutoff", KeyType::Required, cutoff_);

  FANCY_ASSERT(x_range.size() == 2, "Invalid x-range given for Simple Rectangular PV.");
  FANCY_ASSERT(y_range.size() == 2, "Invalid y-range given for Simple Rectangular PV.");
  FANCY_ASSERT(z_range.size() == 2, "Invalid z-range given for Simple Rectangular PV.");
  for(int i = 0; i < 2; i++){
    axis_ranges_[i*3] = x_range[i];
    axis_ranges_[i*3+1] = y_range[i];
    axis_ranges_[i*3+2] = z_range[i];
  }
  return;
}

PV_INDUSRect::PV_INDUSRect(std::array<double,2> xrange, std::array<double,2> yrange, std::array<double,2> zrange)
{
  for(int i = 0; i < 2; i++){
    axis_ranges_[i*3] = xrange[i];
    axis_ranges_[i*3+1] = yrange[i];
    axis_ranges_[i*3+2] = zrange[i];
  }
  return;
}

double PV_INDUSRect::compute(Vec3<double> position){
  Vec3<double> phi;
  for(int i = 0; i < 3; i++){
    phi[i] = h_x(position[i], axis_ranges_[i], axis_ranges_[i+3], sigma_, cutoff_);
  }
  return phi[0]*phi[1]*phi[2];
}