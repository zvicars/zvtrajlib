#include "Calc_IsosurfaceMultiphase.hpp"
#include "../../tools/cellgrid.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/Matrix.hpp"
#include "../../tools/smearfunctions.hpp"
#include "../probevolumes/PV_Trivial.hpp"

RGB convertHexString2RGB(std::string hex_string){
  RGB retVal;
  for(char& c : hex_string){
    c = std::tolower(c);
  }
  std::vector<int> char2int;
  if(hex_string.length() == 8){
    FANCY_ASSERT(hex_string.at(0) == '0' && hex_string.at(1) == 'x', 
                "Hex codes should be 6 characters, 0-9,A-F and can optionally be preceded by a 0x, case insensitive.");
    hex_string.erase(hex_string.begin(), hex_string.begin()+2);
  }
  FANCY_ASSERT(hex_string.length() == 6, "invalid number of characters in hex code string");
  std::array<int, 6> rgb_vals;
  unsigned int counter = 0;
  for(char c : hex_string){
    //check to see if it's 0-9 or a-f
    FANCY_ASSERT(  (c >= (char)48 && c <= (char)57) || (c>= 97 && c <= (char)102), "invalid hex character provided: " + c);
    //case 1, number specified
    if(c >= 48 && c <= 57){
      rgb_vals[counter] = (int)c - 48;
    }
    else rgb_vals[counter] = (int)c - 87;
    counter++;
  }
  retVal.r = rgb_vals[0]*16 + rgb_vals[1];
  retVal.g = rgb_vals[2]*16 + rgb_vals[3];
  retVal.b = rgb_vals[4]*16 + rgb_vals[5];
  //std::cout << hex_string << "  " << retVal.r << "  " << retVal.g << "  " << retVal.b << std::endl;
  //std::cin.get();
  return retVal;
}

Calc_IsosurfaceMultiphase::Calc_IsosurfaceMultiphase(InputPack& input):Calculation{input}
{
  frame_counter_ = 0;
  std::vector<double> npoints;
  input.params().readVector("npoints", KeyType::Required, npoints);
  FANCY_ASSERT(npoints.size() == 3, "Invalid number of dimensions in npoints of isosurface calculation.");
  for(int i = 0; i < 3; i++){
    npoints_[i] = npoints[i];
  }

  std::vector<std::string> agnames;
  input.params().readVector("atom_groups", KeyType::Required, agnames);
  for(const auto& name : agnames){
    AtomGroup* agptr = input.findAtomGroup(name);
    FANCY_ASSERT(agptr != 0, "Failed to find atom group " + name + ".");
    atomgroups_.push_back(agptr);
  }
  agnames_ = agnames;
  num_groups_ = atomgroups_.size();

  std::string pv_name="";
  pv_specified_ = input.params().readString("probe_volume", KeyType::Optional, pv_name);
  if(pv_specified_){
    pv_ = input.findProbeVolume(pv_name);
    FANCY_ASSERT(pv_ != 0, "Failed to find probevolume " + pv_name );
  }
  else{
    pv_ = new PV_Trivial();
  }

  input.params().readVector("distance_rmax", KeyType::Required, distance_rmax_);
  for(auto distance : distance_rmax_){
    FANCY_ASSERT(distance > 0, "Invalid rmax given for isosurface calculation.");
  }

  input.params().readVector("distance_sigmas", KeyType::Required, distance_sigmas_);
  for(auto sigma : sigmas_){
    FANCY_ASSERT(sigma > 0, "Invalid sigma given for isosurface calculation.");
  } 

  input.params().readVector("num_neighbor_threshold", KeyType::Required, num_neighbor_thresholds_);
  for(auto distance : num_neighbor_thresholds_){
    FANCY_ASSERT(distance > 0, "Invalid num_neighbor_threshold given for isosurface calculation.");
  }

  input.params().readVector("num_neighbor_sigmas", KeyType::Required, num_neighbor_sigmas_);
  for(auto sigma : num_neighbor_sigmas_){
    FANCY_ASSERT(sigma > 0, "Invalid sigma given for isosurface calculation.");
  }   

  outputPLY_ = input.params().readVector("hex_colors", KeyType::Optional, hex_colors_);

  input.params().readVector("sigmas", KeyType::Required, sigmas_);
  for(auto sigma : sigmas_){
    FANCY_ASSERT(sigma > 0, "Invalid sigma given for isosurface calculation.");
  }

  input.params().readVector("densities", KeyType::Required, densities_);
  for(auto density : densities_){
     FANCY_ASSERT(density > 0, "Invalid density given for isosurface calculation.");
  }

  input.params().readVector("isovalues", KeyType::Required, isovalues_);
  for(auto isovalue : isovalues_){
    FANCY_ASSERT(isovalue > 0, "Invalid isovalue given for isosurface calculation.");
  }

  method_ = "golosio";
  input.params().readString("method", KeyType::Optional, method_);
  FANCY_ASSERT(method_ == "golosio", "Invalid method chosen for instantaneous interface calculation, valid options are \'golosio\' and \'rchandra\'.");

  return;
}
void Calc_IsosurfaceMultiphase::update(){
  if(hasUpdated()) return;
  Calculation::update();
  for(int i = 0; i < 3; i++){
    box_size_[i] = box->boxvec[i][i];
  }
  if(!initialized_){
    frames_.resize(num_groups_);
    filepaths_.resize(num_groups_);
    meshes_.resize(num_groups_);
    for(int i=0; i<num_groups_; i++){
      //initialize voxelgrids
      frames_[i].initialize(npoints_, box_size_, densities_[i], sigmas_[i], isovalues_[i]);
      //set filepaths for each stla output and delete any existing files with the same name
      filepaths_[i] = base_ + "_" + agnames_[i] + ".stla";
      std::ofstream ofile(filepaths_[i]);
      FANCY_ASSERT(ofile.is_open(), "Failed to open filestream for " + filepaths_[i]);
      ofile.close();
    } 
    initialized_ = 1;
  }
  else{
    for(int i=0; i<num_groups_; i++){
      frames_[i].setLength(box_size_);
      frames_[i].clear();
    } 
  }
  return;
}
void Calc_IsosurfaceMultiphase::calculate(){
  if(!doCalculate()) return;
  //it would be very nice to parallize this, but the add_gaussian() function could lead to data races
  //instead, I placed omp loops inside of add-gaussian function
  int counter = 0;
  for(int i = 0; i < num_groups_; i++){
    for(int j = 0; j < atomgroups_[i]->getIndices().size(); j++ ){
      int idx = atomgroups_[i]->getIndices()[j];
      if(pv_->compute(box->atoms[idx].x) < 0.5) continue;
      frames_[i].add_gaussian(box->atoms[idx].x);
    }
    counter++;
  }
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

void vertexNeighborsRef(const Mesh& mesh, const Mesh& ref_mesh, const CellGrid& ref_cg, 
                        const std::array<double,3>& box_size, const std::array<int,3>& npoints, 
                        double distance_rmax, double distance_sigma, double num_neighbor_cutoff, 
                        double num_neighbor_sigma, std::vector<double>& vertexNeighborCount){
  vertexNeighborCount.resize(mesh.nvtx,0);
  Vec3<double> box_size_offset;
  for(int i=0; i<3; i++){
    box_size_offset[i] = box_size[i] / (double)npoints[i];
  }
  for(std::size_t i = 0; i < (std::size_t)mesh.nvtx; i++){
    const auto& vertex = mesh.vertices[i];
    auto indices = ref_cg.getNearbyIndices(vertex);
    int n_indices = indices.size();
    double stepEval = 0.0;
    for(int j = 0; j < indices.size(); j++){
      Vec3<double> ref_vertex = ref_mesh.vertices[indices[j]];
      getNearestImage3D(ref_vertex, vertex, box_size);
      double distance = norm2(vertex-ref_vertex);
      stepEval += h_r(distance, distance_rmax, distance_sigma, 2.0*distance_sigma);
    }
    vertexNeighborCount[i] += 1.0-h_r(stepEval, num_neighbor_cutoff, num_neighbor_sigma, 2.0*distance_sigma);
  }
  return;
}

void Calc_IsosurfaceMultiphase::output(){
  if(doOutput()){
    //create a buffered version of the voxel grid
    Vec3<int> npoints_temp;
    Vec3<double> box_size_temp;
    for(int i = 0; i < 3; i++){
      npoints_temp[i] = npoints_[i] + 1;
      box_size_temp[i] = box_size_[i] + (box_size_[i] / (double)npoints_[i]);
    }
    for(int i = 0; i < num_groups_; i++){
      VoxelGrid frame_temp(npoints_temp, box_size_temp, densities_[i], sigmas_[i], isovalues_[i]);
      for(int ix = 0; ix < npoints_temp[0]; ix++){
        for(int iy = 0; iy < npoints_temp[1]; iy++){
          for(int iz = 0; iz < npoints_temp[2]; iz++){ 
            frame_temp.setGridVal(ix, iy, iz, frames_[i].getGridVal(wrapIndex(ix, npoints_[0]), wrapIndex(iy, npoints_[1]), wrapIndex(iz, npoints_[2])) );
          }
        }
      }
      //add a buffer to account for periodic boundaries
      marchingCubes(method_, frame_temp, meshes_[i]);
    }
    //using each group as a reference to compute per-group vertex weights 
    std::size_t ng = num_groups_;
    Matrix<double, 2> areaMatrix;
    Matrix<std::vector<double>, 2> vertexStateMatrix;
    areaMatrix.initialize({ng,ng});
    areaMatrix.fill(0);
    vertexStateMatrix.initialize({ng,ng});

    for(std::size_t i = 0; i < ng; i++){
      //create cell grid for group
      CellGrid c_ref(distance_rmax_[i] + 2.0*distance_sigmas_[i], box_size_);
      for(std::size_t j = 0; j < (std::size_t)meshes_[i].nvtx; j++){
        c_ref.addIndexToGrid(j, meshes_[i].vertices[j]);
      }
      for(std::size_t j = 0; j < ng; j++){
        std::array<std::size_t, 2> idx2d = {i,j};
        if(i == j){
          vertexStateMatrix.at(idx2d) = std::vector<double>(meshes_[j].nvtx, 1.0);
          continue;
        }
        std::vector<double> vertexNumNeighbors;
        vertexNeighborsRef(meshes_[j], meshes_[i], c_ref, box_size_, npoints_, distance_rmax_[i], distance_sigmas_[i], 
        num_neighbor_thresholds_[i], num_neighbor_sigmas_[i], vertexNumNeighbors);
        vertexStateMatrix.at(idx2d) = vertexNumNeighbors;
      }
    }
    //compute areas
    for(std::size_t i = 0; i < ng; i++){
      std::array<std::size_t, 2> ii = {i,i};
      for(int j = 0; j < meshes_[i].ntri; j++){
        const auto& triangle = meshes_[i].triangles[j];
        double area = computeTriangleArea(triangle, meshes_[i].vertices);
        for(std::size_t k = 0; k < ng; k++){
          if(i == k) continue;
          std::array<std::size_t, 2> ki = {k,i};
          double weighted_area = 0.0;
          for(int vrt = 0; vrt < 3; vrt++){
            int index = triangle.indices[vrt];
            weighted_area += vertexStateMatrix.at(ki)[index]/3.0;
          }
          weighted_area *= area;
          areaMatrix.at(ki) += weighted_area;
        }
        areaMatrix.at(ii) += area;
      }      
    }
    if(outputPLY_){
      std::vector<RGB> rgb_colors(ng);
      for(int i = 0; i < ng; i++){
        rgb_colors[i]=convertHexString2RGB(hex_colors_[i]);
      }
      for(std::size_t i = 0; i < ng; i++){
        const Mesh& mesh = meshes_[i];
        std::vector< std::vector<double> > vertex_data(ng);
        std::vector<RGB> finalRGBColors(mesh.nvtx);
        for(std::size_t j = 0; j < ng; j++){
          std::array<std::size_t, 2> index = {j,i};
          vertex_data[j] = vertexStateMatrix.read_at(index);
        }
        //vertex data needs to be collapsed and a consensus color picked
        for(int vertex_index = 0; vertex_index < mesh.nvtx; vertex_index++){
          std::vector<std::array<double,3> > final_rgb;
          RGB rgb_average; rgb_average.r = 0.0; rgb_average.g = 0.0; rgb_average.b = 0.0;
          double density_average = 1E-10;
          for(int j = 0; j < ng; j++){
            rgb_average = rgb_average.add(rgb_colors[j].mult(vertex_data[j][vertex_index]));
            density_average += vertex_data[j][vertex_index];
          }
          rgb_average = rgb_average.mult(1.0/density_average);
          finalRGBColors[vertex_index] = rgb_average;
          //std::cout << rgb_average.r << "  " << rgb_average.g << "  " << rgb_average.b << std::endl;
          //std::cin.get();
        }
        std::string ply_output = printPLYWithRGB(mesh, finalRGBColors);
        std::ofstream ofile2(agnames_[i] + "_" + std::to_string((int)box->time) + ".ply");
        ofile2 << ply_output << std::endl;
        ofile2.close();
      }
    }

    //convert matrix to vector
    std::vector<double> output_vec;
    output_vec.reserve(ng*ng); 
    //order is 1,1 1,2 1,3 ... 2,1 2,2 2,3 ..
    //moreover, order is reference mesh then test mesh, so {1,2} will tell you the area of isosurface 2 near isosurface 1
    for(std::size_t i = 0; i < ng; i++){
      for(std::size_t j = 0; j < ng; j++){
        std::array<std::size_t, 2> ij = {i,j};
        output_vec.push_back(areaMatrix.at(ij));
      }
    }

    times_.push_back(box->time);
    avecs_.push_back(output_vec);

    printOutput();
  }  
}

std::string Calc_IsosurfaceMultiphase::printConsoleReport(){
  return "";
}

void Calc_IsosurfaceMultiphase::printOutput(){
  for(int i = 0; i < num_groups_; i++){
    std::string filepath = filepaths_[i];
    std::ofstream ofile(filepath, std::ofstream::app);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface step calculation.");
    std::string output;
    printSTL(meshes_[i], output);
    ofile << output;
    ofile.close();
  }
  return;
};

void Calc_IsosurfaceMultiphase::finalOutput(){
  if(output_freq_ <= 0) return;
  std::ofstream output_ts_data(base_ + "_timeseries.txt");
  FANCY_ASSERT(output_ts_data.is_open(), "Failed to open output stream in Calc_IsosurfaceMultiphase");
  output_ts_data << "#time "; 
  for(int i = 0; i < num_groups_; i++){
    for(int j = 0; j < num_groups_; j++){
      output_ts_data << "[ref " << i << ", test" << j << "]   ";
    }
  }
  output_ts_data << "\n";
  for(int i = 0; i < times_.size(); i++){
    output_ts_data << times_[i] << "  ";
    for(int j = 0; j < avecs_[i].size(); j++){
      output_ts_data << avecs_[i][j] << "  ";
    }
    output_ts_data << "\n";
  }
  output_ts_data.close();
  return;
};