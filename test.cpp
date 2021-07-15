#include "interface/interface.h"
#include <iostream>
int main()
{
  std::cout << "Starting test..." << std::endl;
  XDRTrajectory traj("test.xtc");
  Box b1;
  while(traj.nextFrame()){
    traj.getFrame(b1);
    double pvparms[6] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    int count = OP_Nv(pvparms, b1);
    std::cout << b1.time << "  " << count << std::endl;
  }
  return 0;
}