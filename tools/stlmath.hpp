#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <stdexcept>
template <class T, std::size_t N>
static inline std::array<T, N> operator+(const std::array<T, N>& a, const std::array<T, N>& b){
  std::array<T,N> retVal;
  for(std::size_t i = 0; i < N; i++){
    retVal[i] = a[i] + b[i];
  }
  return retVal;
}

template <class T, std::size_t N>
static inline std::array<T, N> operator-(const std::array<T, N>& a, const std::array<T, N>& b){
  std::array<T,N> retVal;
  for(std::size_t i = 0; i < N; i++){
    retVal[i] = a[i] - b[i];
  }
  return retVal;
}

template <class T, std::size_t N>
static inline std::array<T, N> operator*(const std::array<T, N>& a, const std::array<T, N>& b){
  std::array<T,N> retVal;
  for(std::size_t i = 0; i < N; i++){
    retVal[i] = a[i] * b[i];
  }
  return retVal;
}

template <class T, std::size_t N>
static inline std::array<T, N> operator*(const std::array<T, N>& a, T b){
  std::array<T,N> retVal;
  for(std::size_t i = 0; i < N; i++){
    retVal[i] = a[i] * b;
  }
  return retVal;
}

template <class T, std::size_t N>
static inline std::array<T, N> operator*(T b, const std::array<T, N>& a){
  std::array<T,N> retVal;
  for(std::size_t i = 0; i < N; i++){
    retVal[i] = a[i] * b;
  }
  return retVal;
}

template <class T, std::size_t N>
static inline std::array<T, N> operator/(const std::array<T, N>& a, T b){
  std::array<T,N> retVal;
  for(std::size_t i = 0; i < N; i++){
    retVal[i] = a[i] * b;
  }
  return retVal;
}

template <class T, std::size_t N>
static inline double norm2(const std::array<T, N>& a){
  double retVal = 0.0;
  for(std::size_t i = 0; i < N; i++){
    retVal += a[i]*a[i];
  }
  return std::sqrt(retVal);
}

template <class T>
static inline std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b){
  std::size_t na = a.size();
  std::size_t nb = b.size();
  if(na != nb) throw 7;
  std::vector<T> retVal(na);
  for(std::size_t i = 0; i < na; i++){
    retVal[i] = a[i] + b[i];
  }
  return retVal;
}
template <class T>
static inline std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
  std::size_t na = a.size();
  std::size_t nb = b.size();
  if(na != nb) throw 7;
  std::vector<T> retVal(na);
  for(std::size_t i = 0; i < na; i++){
    retVal[i] = a[i] - b[i];
  }
  return retVal;
}

template <class T>
static inline std::vector<T> operator*(const std::vector<T>& a, T b){
  std::vector<T> retVal(a.size());
  for(std::size_t i = 0; i < a.size(); i++){
    retVal[i] = a[i] * b;
  }
  return retVal;
}

template <class T>
static inline std::vector<T> operator*(T b, const std::vector<T>& a){
  std::vector<T> retVal(a.size());
  for(std::size_t i = 0; i < a.size(); i++){
    retVal[i] = a[i] * b;
  }
  return retVal;
}

template <class T>
static inline double norm2(const std::vector<T>& a){
  double retVal = 0.0;
  std::size_t N = a.size();
  for(std::size_t i = 0; i < N; i++){
    retVal += a[i]*a[i];
  }
  return std::sqrt(retVal);
}

template <class T>
static inline T mean(const std::vector<T>& a){
  double retVal = 0.0;
  for(auto val : a){
    retVal += val;
  }
  retVal /= a.size();
  return retVal;
}

template <class T>
static inline T var(const std::vector<T>& a, T mean){
  double retVal = 0.0;
  for(auto val : a){
    retVal += (val - mean)*(val-mean);
  }
  retVal /= a.size();
  return retVal;
}

template <class T> 
static inline std::array<T,3> cross(const std::array<T,3>& a, const std::array<T,3>& b){
  std::array<T,3> out;
  out[0] = a[1]*b[2] - a[2]*b[1];
  out[1] = a[2]*b[0] - a[0]*b[2];
  out[2] = a[0]*b[1] - a[1]*b[0];
  return out;
}
template <class T> 
static inline T dot(const std::array<T,3>& a, const std::array<T,3>& b){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

template <class T, std::size_t N>
static inline std::array<T, N> vec2Array(std::vector<T> v){
  if(v.size() != N){
    throw std::bad_cast();
  }
  std::array<T,N> ret;
  for(int i = 0; i < N; i++){
    ret[i] = v[i]; 
  }
  return ret;
}