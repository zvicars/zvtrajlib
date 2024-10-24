#pragma once
#include "../interface/datatypes.hpp"
#include "GroManipData.hpp"
#include <cmath>
#include <limits>
class Volume{
  public:
    Volume(){
      return;
    }
    virtual bool isInside(const Vec3<double>& x) const = 0;
    virtual bool isInside(const Vec3<double>& x, const Vec3<double>& boxvec) const {
      //check for each periodic image
      for(int ix = -1; ix <= 1; ix++){
      for(int iy = -1; iy <= 1; iy++){
      for(int iz = -1; iz <= 1; iz++){
            Vec3<double> temp_pos = x;
            temp_pos[0] += ix*boxvec[0]; temp_pos[1] += iy*boxvec[1]; temp_pos[2] += iz*boxvec[2];
            if(isInside(temp_pos)) return 1;
      }}}
      return 0;
    }
    virtual std::vector<double> getBoundingBox() = 0;

};

class UnionVolume : public Volume{
  public:
    UnionVolume(std::vector<std::string> args){
      return;
    };
    UnionVolume(GroManipData& data, std::vector<std::string> args){
      for(auto& name : args){
        auto val = data.findVolume(name);
        FANCY_ASSERT(val != 0, "Failed to find volume in UnionVolume definition.");
        volumes.push_back(val);
      }
      return;
    }
    virtual bool isInside(const Vec3<double>& x) const{
      for(auto& volume : volumes){
        if(volume->isInside(x)) return 1;
      }
      return 0;
    }
    virtual std::vector<double> getBoundingBox(){
      std::vector<double> vec(6, 0.0);
      for(int i = 0; i < 3; i++){
        double xmin = std::numeric_limits<double>::min();
        double xmax = std::numeric_limits<double>::max();
        for(auto vol : volumes){
          auto bb = vol->getBoundingBox();
          xmin = std::min(bb[i], xmin);
          xmax = std::max(bb[i+3], xmax);
        }
        vec[i] = xmin;
        vec[i+3] = xmax;
      }
      return vec;
    }
  protected:
    std::vector<Volume*> volumes;
};

class SphereVolume : public Volume{
  public:
    SphereVolume(std::vector<std::string> args){
      FANCY_ASSERT(args.size() == 4, "Incorrect number of arguments in SphereVolume definition");
      for(int i = 0; i < 3; i++){
        center[i] = std::stod(args[i]);
      }
      radius2 = std::stod(args[3])*std::stod(args[3]);
      return;
    }
    virtual bool isInside(const Vec3<double>& x) const{
      double dist2 = 0.0;
      for(int i = 0; i < 3; i++){
        dist2 +=  (x[i] - center[i]) * (x[i] - center[i]);
      }
      if(dist2 < radius2) return 1;
      return 0;
    }
    virtual std::vector<double> getBoundingBox(){
      std::vector<double> vec(6, 0.0);
      for(int i = 0; i < 3; i++){
        vec[i] = center[i] - std::sqrt(radius2);
        vec[i+3] = center[i] + std::sqrt(radius2);
      }
      return vec;
    }
  protected:
  Vec3<double> center;
  double radius2;
};

class CuboidVolume : public Volume{
  public:
    CuboidVolume(std::vector<std::string> args){
      FANCY_ASSERT(args.size() == 6, "Incorrect number of arguments in CuboidVolume definition");
      for(int i = 0; i < 3; i++){
        min_x[i] = std::stod(args[i]);
        max_x[i] = std::stod(args[i+3]);
      }
      return;
    }
    virtual bool isInside(const Vec3<double>& x) const{
      for(int i = 0; i < 3; i++){
        if(x[i] < min_x[i] || x[i] > max_x[i]) return 0;
      }
      return 1;
    }
    virtual std::vector<double> getBoundingBox(){
      std::vector<double> vec(6, 0.0);
      for(int i = 0; i < 3; i++){
        vec[i] = min_x[i];
        vec[i+3] = max_x[i];
      }
      return vec;
    }
  protected:
    Vec3<double> min_x;
    Vec3<double> max_x;
};

class CylinderVolume : public Volume{
  public:
    CylinderVolume(std::vector<std::string> args){
      FANCY_ASSERT(args.size() == 6, "Incorrect number of arguments in CylinderVolume definition");
      for(int i = 0; i < 3; i++){
        origin[i] = std::stod(args[i]);
      }
      radius2 = std::stod(args[3])*std::stod(args[3]);
      height = std::stod(args[4]);
      axis = std::stoi(args[5]);
      return;
    }
    virtual bool isInside(const Vec3<double>& x) const{
      double dist2 = 0.0;
      for(int i = 0; i < 3; i++){
        if(i == axis) continue;
        dist2 +=  (x[i] - origin[i]) * (x[i] - origin[i]);
      }
      if(dist2 <= radius2){
        if(x[axis] >= origin[axis] && x[axis] <= origin[axis] + height){
          return 1;
        }
      }
      return 0;
    }
    virtual std::vector<double> getBoundingBox(){
      std::vector<double> vec(6, 0.0);
      for(int i = 0; i < 3; i++){
        vec[i] = origin[i] - std::sqrt(radius2);
        vec[i+3] = origin[i] + std::sqrt(radius2);
      }
      vec[axis] = origin[axis];
      vec[axis+3] = origin[axis] + height;
      return vec;
    }
  protected:
  Vec3<double> origin;
  double height;
  int axis;
  double radius2;
};


class EllipseVolume : public Volume{
  public:
    EllipseVolume(std::vector<std::string> args){
      FANCY_ASSERT(args.size() == 6, "Incorrect number of arguments in CylinderVolume definition");
      for(int i = 0; i < 3; i++){
        origin[i] = std::stod(args[i]);
      }
      semi_axes = {std::stod(args[3])*std::stod(args[3]), std::stod(args[4])*std::stod(args[4]), std::stod(args[5])*std::stod(args[5])};
      for(int i = 0; i < 3; i++){
        inv_axes[i] = 1.0/semi_axes[i];
      }
      return;
    }
    virtual bool isInside(const Vec3<double>& x) const{
      double dist2 = 0.0;
      //evaluate ellipse function x2/a2 + y2/b2 ...
      double eval = 0.0;

      for(int i = 0; i < 3; i++){
        eval += std::pow(x[i] - origin[i],2)*inv_axes[i];
      }

      if(eval <= 1.0){
          return 1;
      }

      return 0;
    }
    virtual std::vector<double> getBoundingBox(){
      std::vector<double> vec(6, 0.0);
      for(int i = 0; i < 3; i++){
        vec[i] = origin[i] - std::sqrt(semi_axes[i]);
        vec[i+3] = origin[i] + std::sqrt(semi_axes[i]);
      }
      return vec;
    }
  protected:
  Vec3<double> origin, semi_axes, inv_axes;
};

static inline Volume* createPrimitiveVolume(std::string name, std::vector<std::string> args){
  if(name == "cuboid") return new CuboidVolume(args);
  if(name == "sphere") return new SphereVolume(args);
  if(name == "cylinder") return new CylinderVolume(args);
  if(name == "ellipse") return new EllipseVolume(args);
  return 0;
}

static inline Volume* createUnionVolume(GroManipData& data, std::vector<std::string> args){
  return new UnionVolume(data, args);
  return 0;
}