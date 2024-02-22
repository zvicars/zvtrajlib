#include "Calc_Isosurface3phase.hpp"
#include "../../tools/cellgrid.hpp"
#include "../../tools/stlmath.hpp"

Calc_Isosurface3phase::Calc_Isosurface3phase(InputPack& input):Calculation{input}
{
  frame_counter_ = 0;
  std::vector<double> npoints;
  input.params().readVector("npoints", KeyType::Required, npoints);
  FANCY_ASSERT(npoints.size() == 3, "Invalid number of dimensions in npoints of isosurface calculation.");
  for(int i = 0; i < 3; i++){
    npoints_[i] = npoints[i];
  }

  std::string s_agname, l_agname;
  input.params().readString("solid_atom_group", KeyType::Required, s_agname);
  input.params().readString("liquid_atom_group", KeyType::Required, l_agname);
  solid_atom_group_ = input.findAtomGroup(s_agname);
  liquid_atom_group_ = input.findAtomGroup(l_agname);
  FANCY_ASSERT(solid_atom_group_ != 0, "Failed to find solid atom group.");
  FANCY_ASSERT(liquid_atom_group_ != 0, "Failed to find liquid atom group.");
  
  std::string pv_name;
  input.params().readString("probe_volume", KeyType::Required, pv_name);
  pv_ = input.findProbeVolume(pv_name);
  FANCY_ASSERT(pv_ != 0, "Failed to find probevolume.");

  input.params().readNumber("distance_threshold", KeyType::Required, distance_threshold_);

  input.params().readNumber("liquid_sigma", KeyType::Required, liquid_sigma_);
  input.params().readNumber("solid_sigma", KeyType::Required, solid_sigma_);
  FANCY_ASSERT(liquid_sigma_ > 0, "Invalid sigma given for isosurface calculation.");
  FANCY_ASSERT(liquid_sigma_ > 0, "Invalid sigma given for isosurface calculation.");

  input.params().readNumber("liquid_density", KeyType::Required, liquid_density_);
  FANCY_ASSERT(liquid_density_ > 0, "Invalid density given for isosurface calculation.");
  input.params().readNumber("solid_density", KeyType::Required, solid_density_);
  FANCY_ASSERT(solid_density_ > 0, "Invalid density given for isosurface calculation.");

  input.params().readNumber("liquid_isovalue", KeyType::Required, liquid_isovalue_);
  FANCY_ASSERT(liquid_isovalue_ > 0, "Invalid isovalue given for isosurface calculation.");
  input.params().readNumber("solid_isovalue", KeyType::Required, solid_isovalue_);
  FANCY_ASSERT(solid_isovalue_ > 0, "Invalid isovalue given for isosurface calculation.");


  method_ = "golosio";
  input.params().readString("method", KeyType::Optional, method_);
  FANCY_ASSERT(method_ == "golosio", "Invalid method chosen for instantaneous interface calculation, valid options are \'golosio\' and \'rchandra\'.");
  return;
}
void Calc_Isosurface3phase::update(){
  if(hasUpdated()) return;
  Calculation::update();
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  if(!initialized_){
    liquid_average_.initialize(npoints_, box_size_, liquid_density_, liquid_sigma_, liquid_isovalue_);
    liquid_frame_.initialize(npoints_, box_size_, liquid_density_, liquid_sigma_, liquid_isovalue_);
    solid_average_.initialize(npoints_, box_size_, solid_density_, solid_sigma_, solid_isovalue_);
    solid_frame_.initialize(npoints_, box_size_, solid_density_, solid_sigma_, solid_isovalue_);    

    std::string filepath1 = base_ + "_solidinterface"+ ".stla";
    std::string filepath2 = base_ + "_liquidinterface"+ ".stla";


    initialized_ = 1;
  }
  else{
    liquid_frame_.setLength(box_size_);
    liquid_frame_.clear();
    solid_frame_.setLength(box_size_);
    solid_frame_.clear();
  }
  return;
}
void Calc_Isosurface3phase::calculate(){
  if(!doCalculate()) return;
  //it would be very nice to parallize this, but the add_gaussian() function could lead to data races
  //instead, I placed omp loops inside of add-gaussian function
  for(int i = 0; i < liquid_atom_group_->getIndices().size(); i++ ){
    int idx = liquid_atom_group_->getIndices()[i];
    if(pv_->compute(box->atoms[idx].x) < 0.5) continue;
    liquid_frame_.add_gaussian(box->atoms[idx].x);
  }
  liquid_average_.sumInPlace(liquid_frame_);

  for(int i = 0; i < solid_atom_group_->getIndices().size(); i++ ){
    int idx = solid_atom_group_->getIndices()[i];
    solid_frame_.add_gaussian(box->atoms[idx].x);
  }
  solid_average_.sumInPlace(solid_frame_);

  frame_counter_++;
  return;
}

//triangle holds indices for the vertices
static inline double computeTriangleArea(const Triangle& t1, const std::vector<std::array<double,3> >& vertices){
  //define ab
  Vec3<double> ab = vertices[t1.indices[1]] - vertices[t1.indices[0]];
  Vec3<double> ac = vertices[t1.indices[2]] - vertices[t1.indices[0]];
  return norm2(cross(ab,ac))*0.5;
}

void Calc_Isosurface3phase::output(){
  if(doOutput()){
    marchingCubes(method_, solid_frame_, solid_mesh_);
    marchingCubes(method_, liquid_frame_, liquid_mesh_);
    CellGrid c_solid(distance_threshold_, box_size_);
    for(int i = 0; i < solid_mesh_.nvtx; i++){
      c_solid.addIndexToGrid(i, solid_mesh_.vertices[i]);
    }

    //finds all points on the liquid interface that neighbor the solid interface
    //these vertices will be treated as "solid-neighboring"
    std::vector<bool> liq_neighbors_solid_(liquid_mesh_.nvtx, 0);
    for(std::size_t i = 0; i < liquid_mesh_.nvtx; i++){
      const auto& vertex = liquid_mesh_.vertices[i];
      auto indices = c_solid.getNearbyIndices(vertex);
      int n_indices = indices.size();
      for(int j = 0; j < indices.size(); j++){
        Vec3<double> solid_vertex = solid_mesh_.vertices[indices[j]];
        getNearestImage3D(solid_vertex, vertex, box_size_);
        double distance = norm2(vertex-solid_vertex);
        if(distance <= distance_threshold_) liq_neighbors_solid_[i] = 1;
      }

    }

    //compute the total (and weighted) area of the liquid
    liquid_area_ = 0.0;
    solid_area_ = 0.0;
    double weighted_liquid_area_ = 0.0;

    for(int i = 0; i < liquid_mesh_.ntri; i++){
      const auto& triangle = liquid_mesh_.triangles[i];
      double area = computeTriangleArea(triangle, liquid_mesh_.vertices);
      double weighted_area = 0.0;
      for(int j = 0; j < 3; j++){
        int index = triangle.indices[j];
        if(liq_neighbors_solid_[index] == 1){
          weighted_area += 1.0/3.0;
        }
      }
      weighted_area *= area; //total liquid area
      liquid_area_ += area;
      weighted_liquid_area_ += weighted_area;
    }

    for(int i = 0; i < solid_mesh_.ntri; i++){
      const auto& triangle = solid_mesh_.triangles[i];
      double area = computeTriangleArea(triangle, solid_mesh_.vertices);
      solid_area_ += area;
    }
    time_.push_back(box->time);
    sv_area.push_back(solid_area_ - weighted_liquid_area_);
    sl_area.push_back(weighted_liquid_area_);
    lv_area.push_back(liquid_area_ - weighted_liquid_area_);

    printOutput();
  }  
}

std::string Calc_Isosurface3phase::printConsoleReport(){
  return "";
}

void Calc_Isosurface3phase::printOutput(){
  std::string filepath1 = base_ + "_solidinterface"+ ".stla";
  std::string filepath2 = base_ + "_liquidinterface"+ ".stla";

  std::ofstream ofile(filepath1, std::ofstream::app);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface step calculation.");
  std::string output;
  printSTL(solid_mesh_, output);
  ofile << output;
  ofile.close();
  
  ofile.open(filepath2, std::ofstream::app);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface step calculation.");
  output="";
  printSTL(liquid_mesh_, output);
  ofile << output;
  ofile.close(); 
};

void Calc_Isosurface3phase::finalOutput(){
  if(output_freq_ <= 0) return;
  std::ofstream output_ts_data(base_ + "_timeseries.txt");
  FANCY_ASSERT(output_ts_data.is_open(), "Failed to open output stream in Calc_Isosurface3phase");
  output_ts_data << "#time    sl_area    sv_area    lv_area    \n";
  for(int i = 0; i < time_.size(); i++){
    output_ts_data << time_[i] << "  " << sl_area[i] << "  " << sv_area[i] << "  " << lv_area[i] << "\n";
  }
  output_ts_data.close();
  return;
};