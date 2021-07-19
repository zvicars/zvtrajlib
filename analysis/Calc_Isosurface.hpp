#pragma once
#include "Calculation.hpp"
#include "../cubes/MarchingCubesInterface.hpp"
class Calc_Isosurface : public Calculation{
public:
    Calc_Isosurface(InputPack& input);
    ~Calc_Isosurface();
    virtual void calculate(const Box& box);
    virtual std::string printConsoleReport();    

private:
    std::vector<double> areas_;
    double area_, sigma_;
    Mesh mesh_;
    VoxelGrid* frame_;
    VoxelGrid* average_; //averaging only works for fixed box size
};
