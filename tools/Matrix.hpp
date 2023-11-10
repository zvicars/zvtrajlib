#pragma once
#include <vector>
#include <array>
#include <cassert>
#include <iostream>
#include <iterator>
template <class T, std::size_t N> class Matrix{
public:
  Matrix(){
    return;
  }
  Matrix(std::array<std::size_t, N> dims){
    initialize(dims);
    return;
  }
  void initialize(std::array<std::size_t, N> dims){
    std::size_t vector_size = 1;
    for(auto dim : dims){
      vector_size *= dim;
    }
    data_.resize(vector_size);
    size_ = dims;
    return; 
  }
  void reset(std::array<std::size_t, N> dims){
    data_.clear();
    std::size_t vector_size = 1;
    for(auto dim : dims){
      vector_size *= dim;
    }
    data_.resize(vector_size);
    size_ = dims;
    return; 
  }
  T& at(const std::array<std::size_t, N>& index){
    return data_[mapN1(index)];
  }
  const T& read_at(std::array<std::size_t, N> index) const{
    return read_at_1d(mapN1(index));
  }
  T& at(const std::array<int, N>& index){
    std::array<std::size_t, N> new_index;
    for(int i = 0; i < N; i++){
      new_index[i] = index[i];
    }
    return at_1d(mapN1(new_index));
  }
  const T& read_at(const std::array<int, N>& index) const{
    std::array<std::size_t, N> new_index;
    for(int i = 0; i < N; i++){
      new_index[i] = index[i];
    }
    return read_at_1d(mapN1(new_index));
  }  
  T& at_1d(std::size_t index){
    #ifdef DEBUG
      assert(index[i] < data_.size());
    #endif
    return data_[index];
  }
  const T& read_at_1d(std::size_t index) const{
    #ifdef DEBUG
      assert(index[i] < data_.size());
    #endif
    return data_[index];
  }  
  std::array<std::size_t, N> map1N(std::size_t index) const{
    std::array<std::size_t, N> arr;
    std::size_t index_placeholder = index;
    for(std::size_t i = N; i >= 1; i--){
      std::size_t offset = 1;
      for(std::size_t j = i-1; j >= 1; j--){
        offset *= size_[j-1];
      }
      std::size_t n_offset = index_placeholder / offset;
      arr[i-1] = n_offset;
      index_placeholder -= n_offset*offset;
    }
    return arr;
  }
  std::size_t mapN1(std::array<std::size_t, N> index) const{
    #ifdef DEBUG
    for(std::size_t i = 0; i < N; i++){
      assert(index[i] < size_[i]);
    }
    #endif
    std::size_t idx = 0;
    for(std::size_t i = N; i >= 1; i--){
      std::size_t offset = 1;
      for(std::size_t j = i-1; j >= 1; j--){
        offset *= size_[j-1];
      }
      idx += offset*index[i-1];
    }
    return idx;
  }
  void fill(T value){
    for(auto& val : data_){
      val = value;
    }
    return;
  }
  //N dimensional size of array 
  std::array<std::size_t, N> size() const{
    return size_;
  }
  //one-dimensional size of array
  std::size_t size_1d() const{
    return data_.size();
  }
  //remove all of the elements
  void clear(){
    data_.clear();
    for(auto& val : size_){
      val = 0;
    }
    return;
  }
  //add a matrix and change the current array
  void add_inPlace(Matrix<T,N> b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] += b.at_1d(i);
    }
  }
  //add a constant and change the current array
  void add_inPlace(T b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] += b;
    }
  }
  //subtract by a matrix and change the current array
  void sub_inPlace(Matrix<T,N> b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] -= b.at_1d(i);
    }
  }
  //subtract by a constant and change the current array
  void sub_inPlace(T b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] -= b;
    }
  }
  //multiply by a matrix and change the current array
  void mult_inPlace(Matrix<T,N> b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] *= b.at_1d(i);
    }
  }
  //multiply by a constant and change the current array
  void mult_inPlace(T b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] *= b;
    }
  }
  //divide by a matrix and change the current array
  void div_inPlace(Matrix<T,N> b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] /= b.at_1d(i);
    }
  } 
  //divide by a constant and change the current array
  void div_inPlace(T b){
    #pragma omp parallel for
    for(int i = 0; i < data_.size(); i++){
      data_[i] /= b;
    }
  }
  //checks if both arrays are identical
  bool operator==(const Matrix<T,N>& b) const{
    return data_ == b.data_ &&  size_ == b.size_;
  }
  //element-wise addition between two identically-sized matrices 
  Matrix<T,N> operator+(Matrix<T,N>& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    std::size_t numel2 = b.size_1d();
    assert(numel == numel2);
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) + b.at_1d(i);
    }
    return output;
  }
  //element-wise subtraction between two identically-sized matrices 
  Matrix<T,N> operator-(Matrix<T,N>& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    std::size_t numel2 = b.size_1d();
    assert(numel == numel2);
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) - b.at_1d(i);
    }
    return output;
  }
  //element-wise multiplication between two identically-sized matrices 
  Matrix<T,N> operator*(Matrix<T,N>& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    std::size_t numel2 = b.size_1d();
    assert(numel == numel2);
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) * b.at_1d(i);
    }
    return output;
  }
  //element-wise division between two identically-sized matrices 
  Matrix<T,N> operator/(Matrix<T,N>& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    std::size_t numel2 = b.size_1d();
    assert(numel == numel2);
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) / b.at_1d(i);
    }
    return output;
  }
  //element-wise addition
  Matrix<T,N> operator+(const T& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) + b;
    }
    return output;
  }
  //element-wise subtraction
  Matrix<T,N> operator-(const T& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) - b;
    }
    return output;
  } 
  //element-wise multiplication
  Matrix<T,N> operator*(const T& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) * b;
    }
    return output;
  }
  //element-wise division
  Matrix<T,N> operator/(const T& b){
    Matrix<T,N> output(size_);
    std::size_t numel = data_.size();
    #pragma omp parallel for
    for(int i = 0; i < numel; i++){
      output.at_1d(i) = at_1d(i) / b;
    }
    return output;
  }
  //returns a pointer to the start of the vector data
  //use this for 1d sums
  T* data(const T& b){
    return data_.data();
  }
  std::vector<T>& vector(const T& b){
    return data_;
  }  
protected:
  std::vector<T> data_;
  std::array<std::size_t, N> size_;
};


template <class T, std::size_t N, std::size_t M> 
class Matrix2d : public Matrix<T,2>{
public:
  Matrix2d(){
    std::array<std::size_t, 2> dims = {N,M};
    Matrix<T,2>::initialize(dims);
  }
};