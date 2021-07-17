#pragma once
#include "../io/InputParser.hpp"
#include "probevolumes.hpp"
#include "calculations.hpp"
//generic parser type that will return a type T pointer contingent upon the keys and constructor functions registered
template <class T>
class Parser{
public:
  ~Parser{
    for (auto iter = map_.rbegin(); iter != map_.rend(); ++iter) {
        auto ptr = iter.second;
        map.erase(iter);
        delete ptr; 
    }
    return;
  }
  register( std::string key, T* (*f_in)(std::string) ){
    map_.insert({key, f_in});
  };
  T* construct(std::string key, std::string input_args){
    auto it = map_.find(key);
    if(it == map_.end()){
      std::cout << "Key \"" << key << "\" does not have a registered constructing function.";
      throw 0;
    }
    auto f = it->second;
    return f(input_args);
  };
private:
  std::map<std::string, T* (*)(std::string)> map_;
};

//an object registry will be created for each type T, this will map object names to their appropriate T* pointer
//each object registry will have a parser constructed within it, so constructor registration happens via the object register
//all pointers will be deleted by the registry
template <class T>
class ObjectRegistry{
public:
  ~ObjectRegistry{
    for (auto iter = map_.rbegin(); iter != map_.rend(); ++iter) {
        auto ptr = iter.second;
        map.erase(iter);
        delete ptr; 
    }
    return;
  }
  registerParser( std::string key, T* (*f_in)(std::string) ){ //this is how you add derived classes in
    parser_.register(key, f_in);
  }
  registerObject(std::string input_data){
      std::string name, type;
      viaKey<std::string>("name", input_data, name);
      viaKey<std::string>( "type", input_data, type);
      T* obj = parser_.construct(type, input_data);
      map_.insert(name, obj);
  }
  T* findObject(std::string name){
    auto i = map_.find(name);
    if(it == map_.end()){
      std::cout << "Key \"" << key << "\" does not have a registered constructing function.";
      throw 0;
    }
    return it->second;
  }
private:
  Parser<T> parser_;
  std::map<std::string, T*> map_ //a map of all object names to their corresponding base-class pointer
};

//some program state options
namespace reg{
  ObjectRegistry<ProbeVolume> probe_volumes; //as there are many variants of probe volumes and probe volumes may be reused, it's worth registering
  ObjectRegistry<Calculation> calculations; //something that operates on a frame and outputs some data, implementation gets full control on data output
}