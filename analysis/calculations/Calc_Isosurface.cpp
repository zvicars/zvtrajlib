#include "Calc_Isosurface.hpp"
Calc_Isosurface::Calc_Isosurface(InputPack& input):Calculation{input}
{
  frame_counter_ = 0;
  std::vector<double> npoints;
  input.params().readVector("npoints", KeyType::Required, npoints);
  FANCY_ASSERT(npoints.size() == 3, "Invalid number of dimensions in npoints of isosurface calculation.");
  for(int i = 0; i < 3; i++){
    npoints_[i] = npoints[i];
  }

  std::string agname;
  input.params().readString("atom_group", KeyType::Required, agname);
  atom_group_ = input.findAtomGroup(agname);
  FANCY_ASSERT(atom_group_ != 0, "Failed to find specified atom group.");
  
  input.params().readNumber("sigma", KeyType::Required, sigma_);
  FANCY_ASSERT(sigma_ > 0, "Invalid sigma given for isosurface calculation.");
  input.params().readNumber("density", KeyType::Required, density_);
  FANCY_ASSERT(density_ > 0, "Invalid density given for isosurface calculation.");
  input.params().readNumber("isovalue", KeyType::Required, isovalue_);
  FANCY_ASSERT(isovalue_ > 0, "Invalid isovalue given for isosurface calculation.");
  method_ = "golosio";
  input.params().readString("method", KeyType::Optional, method_);
  FANCY_ASSERT(method_ == "golosio", "Invalid method chosen for instantaneous interface calculation, valid options are \'golosio\' and \'rchandra\'.");
  return;
}
void Calc_Isosurface::update(){
  Calculation::update();
  Vec3<double> box_size;
  for(int i = 0; i < 3; i++){
    box_size[i] = box->boxvec[i][i];
  }
  if(!initialized_){
    average_.initialize(npoints_, box_size, density_, sigma_, isovalue_);
    frame_.initialize(npoints_, box_size, density_, sigma_, isovalue_);
    initialized_ = 1;
  }
  else{
    frame_.setLength(box_size);
    frame_.clear();
  }
  return;
}
void Calc_Isosurface::calculate(){
  if(!doCalculate()) return;
  //it would be very nice to parallize this, but the add_gaussian() function could lead to data races
  //instead, I placed omp loops inside of add-gaussian function
  for(int i = 0; i < atom_group_->getIndices().size(); i++ ){
    int idx = atom_group_->getIndices()[i];
    frame_.add_gaussian(box->atoms[idx].x);
  }
  average_.sumInPlace(frame_);
  frame_counter_++;
  return;
}

void Calc_Isosurface::output(){
  if(doOutput()){
    marchingCubes(method_, frame_, mesh_);
    printOutput();
  }  
}

std::string Calc_Isosurface::printConsoleReport(){
  return "";
}

void Calc_Isosurface::printOutput(){
  std::string filepath = base_ + "_interface"+ ".stla";
  std::ofstream ofile(filepath, std::ofstream::app);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface step calculation.");
  std::string output;
  printSTL(mesh_, output);
  ofile << output;
  ofile.close();
};
void Calc_Isosurface::finalOutput(){
  if(output_freq_ <= 0) return;
  average_.scalarMult(1.0/(double)frame_counter_);
  marchingCubes(method_, average_, mesh_);
  std::string filepath = base_ + "_average.ply";
  std::ofstream ofile(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface average.");
  std::string output;
  //printSTL(mesh_, output);
  std::vector<std::array<double, 2> > curv;
  curv = computeMeshCurvature(mesh_, 3);
  std::vector<double> ac(curv.size()), gc(curv.size());
  for(int i = 0; i < curv.size(); i++){
    ac[i] = 0.5*(curv[i][0]+curv[i][1]);
    gc[i] = curv[i][0]*curv[i][1];
  }
  std::string ac_info = printPLYWithCurvature(mesh_, ac); 
  std::string gc_info = printPLYWithCurvature(mesh_, gc);   
  ofile << ac_info;
  ofile.close(); 
  filepath = base_ + "_gaussian.ply";
  ofile.open(filepath);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface average.");
  ofile << gc_info;
  ofile.close();
  return;
};