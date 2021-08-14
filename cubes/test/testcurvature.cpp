#include "../MarchingCubesInterface.hpp"
#include "../meshfiles/objfile.hpp"
#include <string>
#include <fstream>
#include <streambuf>
int main(){
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
  std::vector< std::array<double,2> > curvatures;
  std::string output_info = printPLYWithCurvature(m1, curvatures, 3);
  std::ofstream ofile("TestSphereOutput.ply");
  ofile << output_info;
  ofile.close();
  double sum = 0;
  ofile.open("curvatures.txt");
  for(int i = 0; i < curvatures.size(); i++){
    ofile << curvatures[i][0] << "     " << curvatures[i][1] << std::endl;
    sum += 0.5*(curvatures[i][0] + curvatures[i][1]);
  }
  ofile.close();
  sum *= 1.0/(double)curvatures.size();
  if(fabs(sum + 0.333) < 1e-2) return 0;
  std::cout << "Curvature is not what is expected for example spherical mesh." << std::endl;
  return 1 ;
}