#include "../MarchingCubesInterface.hpp"
#include "../meshfiles/objfile.hpp"
#include <string>
#include <fstream>
#include <streambuf>

bool isNearRel(double actual, double expected, double tol){
  if(actual == 0.0){
    return fabs(expected) < tol;
  }
  return fabs((actual - expected)/actual) < tol;
}

void getCurvatures(const std::vector< std::array<double, 2> >& curvatures, double& av, double& gauss){
  av = 0.0;
  gauss = 0.0;
  for(int i = 0; i < curvatures.size(); i++){
    av += 0.5*(curvatures[i][0] + curvatures[i][1]);
    gauss += curvatures[i][0]*curvatures[i][1];
  }  
  av *= 1.0/(double)curvatures.size();
  gauss *= 1.0/(double)curvatures.size();  
  return;
}

int testSphereCurvature()
{
  std::ifstream t("testfiles/TestSphere.obj");
  if(!t.is_open()){
    std::cout << "Failed to open file." << std::endl;
    return 1;
  }
  std::string str((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());
  t.close();
  ObjFile myObj;
  myObj.load(str);
  Mesh m1 = myObj.getMesh();
  std::vector< std::array<double,2> > curvatures = computeMeshCurvature(m1, 3);
  std::vector<double> ac(m1.nvtx), gc(m1.nvtx);
  for(int i = 0; i < m1.nvtx; i++){
    ac[i] = 0.5*(curvatures[i][0]+curvatures[i][1]);
    gc[i] = curvatures[i][0]*curvatures[i][1];
  }
  std::string output_info = printPLYWithCurvature(m1, ac);
  std::ofstream ofile("TestSphereOutput.ply");
  ofile << output_info;
  ofile.close();
  double av, gauss;
  getCurvatures(curvatures, av, gauss);
  std::cout << "Mean curvature of TestSphere.obj - " << "Observed: " << av << " Expected: " << -0.333 << std::endl;
  std::cout << "Gaussian curvature of TestSphere.obj - " << "Observed: " << gauss << " Expected: " << 1.0/9.0 << std::endl;
  if(isNearRel(-1.0/3.0, av, 0.01) && isNearRel(1.0/9.0, gauss, 0.01)) return 0;
  return 1;
}
int testCylinderCurvature()
{
  std::ifstream t("testfiles/TestCylinder.obj");
  if(!t.is_open()){
    std::cout << "Failed to open file." << std::endl;
    return 1;
  }
  std::string str((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());
  t.close();
  ObjFile myObj;
  myObj.load(str);
  Mesh m1 = myObj.getMesh();
  std::vector< std::array<double,2> > curvatures = computeMeshCurvature(m1, 3);
  std::vector<double> ac(m1.nvtx), gc(m1.nvtx);
  for(int i = 0; i < m1.nvtx; i++){
    ac[i] = 0.5*(curvatures[i][0]+curvatures[i][1]);
    gc[i] = curvatures[i][0]*curvatures[i][1];
  }
  std::string output_info = printPLYWithCurvature(m1, ac);
  std::ofstream ofile("TestCylinderOutput.ply");
  ofile << output_info;
  ofile.close();
  double av, gauss;
  getCurvatures(curvatures, av, gauss);
  std::cout << "Mean curvature of TestCylinder.obj - " << "Observed: " << av << " Expected: " << -1.0/6.0 << std::endl;
  std::cout << "Gaussian curvature of TestCylinder.obj - " << "Observed: " << gauss << " Expected: " << 0.0 << std::endl; 
  if(!isNearRel(-1.0/6.0, av, 0.02)) return 1;
  if(!isNearRel(0.0, gauss, 0.02)) return 1;
  return 0;
}

int testTorusCurvature()
{
  std::ifstream t("testfiles/TestTorus.obj");
  if(!t.is_open()){
    std::cout << "Failed to open file." << std::endl;
    return 1;
  }
  std::string str((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());
  t.close();
  ObjFile myObj;
  myObj.load(str);
  Mesh m1 = myObj.getMesh();
  std::vector< std::array<double,2> > curvatures = computeMeshCurvature(m1, 3);
  std::vector<double> ac(m1.nvtx), gc(m1.nvtx);
  for(int i = 0; i < m1.nvtx; i++){
    ac[i] = 0.5*(curvatures[i][0]+curvatures[i][1]);
    gc[i] = curvatures[i][0]*curvatures[i][1];
  }
  std::string ac_info = printPLYWithCurvature(m1, ac); 
  std::string gc_info = printPLYWithCurvature(m1, gc); 
  std::ofstream ofile("TestTorusOutputAverage.ply");
  ofile << ac_info;
  ofile.close();
  ofile.open("TestTorusOutputGaussian.ply");
  ofile << gc_info;
  ofile.close();
  double av, gauss;
  getCurvatures(curvatures, av, gauss);
  std::cout << "Mean curvature of TestTorus.obj - " << "Observed: " << av << std::endl;
  std::cout << "Gaussian curvature of TestTorus.obj - " << "Observed: " << gauss << std::endl; 
  //if(!isNearRel(-1.0/6.0, av, 0.02)) return 1;
  //if(!isNearRel(0.0, gauss, 0.02)) return 1;
  return 0;
}

int testCubeCurvature()
{
  std::ifstream t("testfiles/TestCube.obj");
  if(!t.is_open()){
    std::cout << "Failed to open file." << std::endl;
    return 1;
  }
  std::string str((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());
  t.close();
  ObjFile myObj;
  myObj.load(str);
  Mesh m1 = myObj.getMesh();
  std::vector< std::array<double,2> > curvatures = computeMeshCurvature(m1, 3);
  std::vector<double> ac(m1.nvtx), gc(m1.nvtx);
  for(int i = 0; i < m1.nvtx; i++){
    ac[i] = 0.5*(curvatures[i][0]+curvatures[i][1]);
    gc[i] = curvatures[i][0]*curvatures[i][1];
  }
  std::string ac_info = printPLYWithCurvature(m1, ac); 
  std::string gc_info = printPLYWithCurvature(m1, gc); 
  std::ofstream ofile("TestCubeOutputAverage.ply");
  ofile << ac_info;
  ofile.close();
  ofile.open("TestCubeOutputGaussian.ply");
  ofile << gc_info;
  ofile.close();
  double av, gauss;
  getCurvatures(curvatures, av, gauss);
  std::cout << "Mean curvature of TestCube.obj - " << "Observed: " << av << std::endl;
  std::cout << "Gaussian curvature of TestCube.obj - " << "Observed: " << gauss << std::endl; 
  //if(!isNearRel(-1.0/6.0, av, 0.02)) return 1;
  //if(!isNearRel(0.0, gauss, 0.02)) return 1;
  return 0;
}

int main(){
  if(testSphereCurvature()) return 1;
  if(testCylinderCurvature()) return 1;
  //test saving feature of obj object
  return 0 ;
}