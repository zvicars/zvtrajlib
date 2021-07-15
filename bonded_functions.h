#pragma once
#include <vector>
#include <array>

template <class T>
using Vec3 = std::array<T, 3>;
template <class T>
using Vec = std::vector<T>;
/*
Pair potentials
*/
//parameters position i, position j, parameters [k, r], force output
inline double harmonic_bond(const Vec3<double>& x1, const Vec3<double>& x2, const Vec<double>& params){
  return 0;
}
//parameters position i, position j, parameters [De, a, re], force output
inline double morse_bond(const Vec3<double>& x1, const Vec3<double>& x2, const Vec<double>& params){
  return 0;
};
//parameters position i, position j, parameters [sigma, epsilon], force output
inline double lj(const Vec3<double>& x1, const Vec3<double>& x2, const Vec<double>& params){
  return 0;
}
/*
Angle potentials
*/
//parameters position i, position j, parameters [k, theta], force output
inline double harmonic_angle(const Vec3<double>& x1, const Vec3<double>& x2, const Vec3<double>& x3, Vec<double>& params){
  return 0;
}

/*
Dihedral potentials
*/
inline double dihedral_angle(const Vec3<double>& x1, const Vec3<double>& x2, const Vec3<double>& x3, Vec<double>& params){
  return 0;
}
