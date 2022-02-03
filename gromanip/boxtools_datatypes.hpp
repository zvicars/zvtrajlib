#pragma  once
#include "../tools/Assert.hpp"
struct SingleInst{
  int i;
  virtual std::string print() = 0;
};

struct SingleType{
  std::string i;
  virtual std::string print() = 0;
};

struct BTAtomType{
  std::string name, type;
  double charge, mass;
  virtual std::string print(){
    std::string output = "";
    output = name + "   " + type + "   " + std::to_string(charge) + "   " + std::to_string(mass);
    return output;
  }
  virtual bool read(std::string in){
    //;name type charge mass
    std::stringstream ss(in);
    ss >> name >> type >> charge >> mass;
    return ss.fail();
  } 
};
struct AtomInst{
  int nr, resnr, cgnr;
  std::string name, type, res;
  double charge, mass;
  AtomInst(int nr_in, int resnr_in, int cgnr_in, const std::string& res_in, const BTAtomType& ref){
    mass = ref.mass;
    charge = ref.charge;
    name = ref.name;
    type = ref.type;
    nr = nr_in;
    res = res_in;
    resnr = resnr_in;
    cgnr = cgnr_in;
    return;
  }
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(nr) + "   " + type + "   " + 
    std::to_string(resnr) + "   " + res + "   " + name + 
    "   " + std::to_string(cgnr) + "   " + std::to_string(charge) + 
    "   " + std::to_string(mass);
    return output;
  }
  virtual bool read(std::string in){
    //;name type charge mass
    std::stringstream ss(in);
    ss >> nr >> type >> resnr >> res >> name >> cgnr >> charge >> mass;
    return ss.fail();
  } 
};

struct PosResType{
  std::string name;
  int function_type;
  double fcx, fcy, fcz;
  virtual std::string print(){
    std::string output = "";
    output = name + "   " + std::to_string(function_type) + "   " + std::to_string(fcx) + "   " + std::to_string(fcy) + "   " + std::to_string(fcz);
    return output;
  }
  virtual bool read(std::string in){
    //;name type charge mass
    std::stringstream ss(in);
    ss >> name >> function_type >> fcx >> fcy >> fcz;
    return ss.fail();
  } 
};
struct PosResInst{
  int i;
  int function_type;
  double fcx, fcy, fcz;
  PosResInst(int a, const PosResType& ref){
    i = a;
    function_type = ref.function_type;
    fcx = ref.fcx;
    fcy = ref.fcy;
    fcz = ref.fcz;
    return;
  }
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(i) + "   " + std::to_string(function_type) + "   " + std::to_string(fcx) + "   " + std::to_string(fcy) + "   " + std::to_string(fcz);
    return output;
  }
  virtual bool read(std::string in){
    //;name type charge mass
    std::stringstream ss(in);
    ss >> i >> function_type >> fcx >> fcy >> fcz;
    return ss.fail();
  } 
};

struct PairType{
  std::string i, j;
  virtual std::string print(){
    std::string output = "";
    output = i + "   " + j;
    return output;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j;
    return ss.fail();
  }
};
struct PairInst{
  int i, j;
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(i) + "   " + std::to_string(j);
    return output;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j;
    return ss.fail();
  }
};

struct BondType : public PairType{
  std::array<double,2> params;
  int function_type;
  virtual std::string print(){
    std::string output = "";
    output = i + "   " + j + "   " + std::to_string(function_type) + "   ";
    for(auto param : params){
      output += std::to_string(param) + "   ";
    }
    return output;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> function_type >> params[0] >> params[1];
    return ss.fail();
  }
};
struct BondInst : public PairInst{
  std::array<double,2> params;
  int function_type;
  BondInst(int a, int b, const BondType& ref){
    i = a; j = b;
    function_type = ref.function_type;
    params = ref.params;
    return;
  }
  BondInst(){
    return;
  }
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(i) + "   " + std::to_string(j) + "   " + std::to_string(function_type) + "   ";
    for(auto param : params){
      output += std::to_string(param) + "   ";
    }
    return output;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> function_type >> params[0] >> params[1];
    return ss.fail();
  }
};

struct PBCBondType : public PairType{
  std::array<double,2> params;
  int function_type;
  Vec3<bool> pbc;
  virtual bool read(std::string in){
    std::stringstream ss(in);
    std::string pbc_string;
    ss >> i >> j >> function_type >> pbc_string >> params[0] >> params[1];
    if(pbc_string.find('x') != std::string::npos) pbc[0] = 1;
    else pbc[0] = 0;
    if(pbc_string.find('y') != std::string::npos) pbc[1] = 1;
    else pbc[1] = 0;
    if(pbc_string.find('z') != std::string::npos) pbc[2] = 1;
    else pbc[2] = 0;    
    return ss.fail();
  }
  virtual std::string print(){
    std::string output = "";
    output = i + "   " + j + "   " + std::to_string(function_type) + "   ";
    for(auto param : params){
      output += std::to_string(param) + "   ";
    }
    return output;
  }
};

struct ExclusionType : public SingleType{
  std::vector<std::string> exclusion_names;
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i;
    std::string entry;
    while(ss >> entry){
      exclusion_names.push_back(entry);
    }
    if(exclusion_names.size() == 0) return 1;
    return 0;
  }
  virtual std::string print(){
    std::string output = "";
    output = i + "   ";
    for(auto name : exclusion_names){
      output += name + "   ";
    }
    return output;
  }
};
struct ExclusionInst : public SingleInst{
  std::vector<int> exclusion_names;
  ExclusionInst(int a, std::vector<int> exclusions){
    i = a;
    exclusion_names = exclusions;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i;
    int entry;
    while(ss >> entry){
      exclusion_names.push_back(entry);
    }
    if(exclusion_names.size() == 0) return 1;
    return 0;
  }
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(i) + "   ";
    for(auto name : exclusion_names){
      output += std::to_string(name) + "   ";
    }
    return output;
  }
};
struct ConstraintType : public PairType{
  int function_type;
  double b0;
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> function_type >> b0;
    return ss.fail();
  }
  virtual std::string print(){
    std::string output = "";
    output = i + "   " + j + "   " + std::to_string(function_type) + "   " + std::to_string(b0);
    return output;
  }
};
struct ConstraintInst : public PairInst{
  int function_type;
  double b0;
  ConstraintInst(int a, int b, const ConstraintType& ref){
    i = a; j = b;
    b0 = ref.b0;
    function_type = ref.function_type;
    return;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> function_type >> b0;
    return ss.fail();
  }
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(i) + "   " + std::to_string(j) + "   " + std::to_string(function_type) + "   " + std::to_string(b0);
    return output;
  }
};
struct AngleType{
  std::string i, j, k;
  std::array<double, 2> params;
  int function_type;
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> k >> function_type >> params[0] >> params[1];
    return ss.fail();
  }
  virtual std::string print(){
    std::string output = "";
    output = i + "   " + j + "   " + k + "   " + std::to_string(function_type) + "   ";
    for(auto param : params){
      output += std::to_string(param) + "   ";
    }
    return output;
  }
};
struct AngleInst{
  int i, j, k;
  std::array<double, 2> params;
  int function_type;
  AngleInst(int a, int b, int c, const AngleType& ref){
    i = a; j = b; k = c;
    function_type = ref.function_type;
    params = ref.params;
    return;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> k >> function_type >> params[0] >> params[1];
    return ss.fail();
  }
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(i) + "   " + std::to_string(j) + "   " + std::to_string(k) + "   " + std::to_string(function_type) + "   ";
    for(auto param : params){
      output += std::to_string(param) + "   ";
    }
    return output;
  }
};

struct DihedralType{
  std::string i, j, k, l;
  std::vector<double> params;
  int function_type;
};





struct Vsite3Type {
  std::string i, j, k, l;
  int function_type;
  std::array<double, 2> params;
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> k >> l >> function_type >> params[0] >> params[1];
    return ss.fail();
  }
  virtual std::string print(){
    std::string output = "";
    output = i + "   " + j + "   " + k + "   " + l + "   " + std::to_string(function_type) + "   ";
    for(auto param : params){
      output += std::to_string(param) + "   ";
    }
    return output;
  }
};
struct Vsite3Inst {
  int i, j, k, l;
  int function_type;
  std::array<double, 2> params;
  Vsite3Inst(int a, int b, int c, int d, const Vsite3Type& ref){
    i = a; j = b; k = c; l = d;
    function_type = ref.function_type;
    params = ref.params;
    return;
  }
  virtual bool read(std::string in){
    std::stringstream ss(in);
    ss >> i >> j >> k >> l >> function_type >> params[0] >> params[1];
    return ss.fail();
  }
  virtual std::string print(){
    std::string output = "";
    output = std::to_string(i) + "   " + std::to_string(j) + "   "
             + std::to_string(k) + "   " + std::to_string(l) + "   " 
             + std::to_string(function_type) + "   ";
    for(auto param : params){
      output += std::to_string(param) + "   ";
    }
    return output;
  }
};