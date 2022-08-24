#include "Matrix.hpp"
#include "stlmath.hpp"

template <class T, std::size_t N>
static Matrix<std::array<T,N>, N> gradient(Matrix<T,N>& input){
  Matrix<std::array<T,N>,N> output(input.size());
  //for each element in the matrix, do a finite-difference estimate in all N directions
  auto numel = input.size_1d();
  for(int i = 0; i < numel; i++){
    //get the N-dimensional index
    std::array<std::size_t,N> index_step = input.map1N(i);
    std::array<T,N> diff_step;
    for(int j = 0; j < N; j++){
      auto index_step_temp = index_step;
      auto forward_step = index_step;
      if(index_step_temp[j] == input.size()[j]-1){
        index_step_temp[j]--;
      }
      else{
        forward_step[j]++;
      }
      diff_step[j] = input.at(forward_step) - input.at(index_step_temp);
    }
    output.at_1d(i) = diff_step;
    //output.at_1d(i) = {1.0,1.0,1.0};
  }
  return output;
}

template <class T, std::size_t N>
static Matrix<Matrix2d<T,N,N>, N> Hessian(Matrix<T,N>& input){
  //each element in the voxel grid will have an associated Hessian
  Matrix<Matrix2d<T,N,N>, N> ret;
  auto size = input.size();
  ret.initialize(size);
  std::array<std::size_t, 2> dim = {N,N};
  //NxN matrix with entries d^2f/dxidxj
  //loop through each element in the array
  #pragma omp parallel for
  for(std::size_t i = 0; i < input.size_1d(); i++){
    Matrix2d<T,N,N> entryTemplate;
    std::array<std::size_t, N> index = input.map1N(i);
    for(std::size_t i2 = 0; i2 < N; i2++){
      for(std::size_t j = 0; j < N; j++){
        std::size_t bounds_i = size[i2], bounds_j = size[j]; 
        auto index_temp = index;
        auto index_shifti = index;
        auto index_shiftj = index;
        auto index_shiftij = index;
        //backward difference at boundaries
        if(index_temp[i2] == bounds_i-1){
          index_temp[i2]--;
        }
        else{
          index_shifti[i2]++;
          index_shiftij[i2]++;
        }
        if(index_temp[j] == bounds_j-1){
          index_temp[j]--;
        }
        else{
          index_shiftj[j]++;
          index_shiftij[j]++;
          //techically wrong, but only on the edges
          if(index_shiftj[j] == bounds_j) index_shiftj[j]--;
          if(index_shiftij[j] == bounds_j) index_shiftij[j]--;
        }
        T diff = input.at(index_shiftij) - input.at(index_shifti) - input.at(index_shiftj) + input.at(index_temp);
        std::array<std::size_t, 2> index2 = {i2,j};
        entryTemplate.at(index2) = diff;
      }
    }
    ret.at_1d(i) = entryTemplate;
  }
  return ret;
}