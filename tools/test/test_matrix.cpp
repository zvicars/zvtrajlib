#include "../Matrix.hpp"
#include "../stlmath.hpp"
#include "../Assert.hpp"
#include <iostream>
int main(){
  //create a 3d matrix that is 5x4x3
  std::cout << "Testing to see if matrix creation, writing, and reading operations are consistent." << std::endl;
  std::array<std::size_t,3> matrix_size = {6,5,4};
  Matrix<double, 3> new_mat(matrix_size);
  std::cout << "Initialized matrix" << std::endl;
  new_mat.fill(0.0);
  std::cout << "Filled matrix" << std::endl;
  std::array<std::size_t,3> index = {1,2,3};
  new_mat.at(index) = 5.0;
  std::cout << "Assigned value" << std::endl;
  for(std::size_t i = 0; i < new_mat.size_1d(); i++){
    auto val = new_mat.at_1d(i);
    if(val == 5.0){
      auto index = new_mat.map1N(i);
      FANCY_ASSERT(new_mat.map1N(i) == index, "Read/write operations are inconsistent. Value written at [1 2 3] was read at [" +
      std::to_string(index[0]) + " " + std::to_string(index[1]) + " " + std::to_string(index[2]) + "]");
    }
  }
  std::cout << "Testing in-place addition, multiplication, addition, and subtraction" << std::endl;
  Matrix<double,3> new_mat2(matrix_size);
  Matrix<double,3> new_mat3(matrix_size);
  Matrix<double,3> new_mat4(matrix_size);
  new_mat4.fill(1.0);
  for(std::size_t i = 0; i < matrix_size[0]; i++){
    for(std::size_t j = 0; j < matrix_size[1]; j++){
      for(std::size_t k = 0; k < matrix_size[2]; k++){
        std::array<std::size_t, 3> index = {i,j,k};
        new_mat2.at(index) = i+j+k;
        new_mat3.at(index) = i+j+k+1;
      }
    }
  }
  FANCY_ASSERT(new_mat4 == (new_mat3 - new_mat2), "Incorrect result for subtraction.");
  FANCY_ASSERT(new_mat2*2.0 == (new_mat2*4.0)*0.5, "Incorrect result for multiplication");

  for(std::size_t i = 0; i < matrix_size[0]; i++){
    for(std::size_t j = 0; j < matrix_size[1]; j++){
      for(std::size_t k = 0; k < matrix_size[2]; k++){
        std::cout << "i j k: " << i << " " << j << " " << k << "\n";
        auto a = new_mat.mapN1({i,j,k});
        auto b = new_mat.map1N(a);
        std::cout << "i j k mapped: " << b[0] << " " << b[1] << " " << b[2] << "\n";
        auto c = new_mat.mapN1(b);
        std::cout << "1d indices " << a << "   " << c << std::endl;
        FANCY_ASSERT(a == c, "coordinate mapping appears to be non-invertible");
  }}}    
  return 0;
}