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
        if(!volume->isInside(x)) return 0;
      }
      return 1;
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

static inline Volume* createPrimitiveVolume(std::string name, std::vector<std::string> args){
  if(name == "cuboid") return new CuboidVolume(args);
  if(name == "sphere") return new SphereVolume(args);
  return 0;
}

static inline Volume* createUnionVolume(GroManipData& data, std::vector<std::string> args){
  return new UnionVolume(data, args);
  return 0;
}