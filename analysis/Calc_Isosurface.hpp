#pragma once
#include "Calculation.hpp"
#include "../cubes/MarchingCubesInterface.hpp"
class Calc_Isosurface : public Calculation{
public:
    Calc_Isosurface(InputPack& input);
    ~Calc_Isosurface(){
        delete frame_;
        delete average_;
    }
    virtual void calculate(const Box& box);
    virtual std::string printConsoleReport();
private:
    int frame_counter_;
    std::vector<double> areas_;
    Vec3<int> npoints_;
    double area_, sigma_, density_, isovalue_;
    Mesh mesh_;
    VoxelGrid* frame_;
    VoxelGrid* average_;
};
