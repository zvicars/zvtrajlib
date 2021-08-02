#include "../VoxelGrid.hpp"
//Testing script for the voxelgrid object type to ensure that it is quantitatively accurate. Two things must remain true, the total density must add up to N/density and the 
//center of mass of identical mass particles must be conserved
//to this end, a single atom will be added to the box in an arbitrary location and these quantities will be computed
//one caveat is that the com might be split across pbcs, so the added particle will be far away from the boundaries, this isn't ideal as an error in pbc treatment would go unnoticed, but
//I want to get this test up and running asap and this covers typical cases
int main(){
    double testing_threshold = 1e-4;
    std::array<int,3> size = {10, 15, 20};
    std::array<double, 3> box_size = {20.0, 15.0, 10.0};
    double density = 1;
    double sigma = 0.5;
    double isovalue = 0.5;
    VoxelGrid v1(size, box_size, density, sigma, isovalue);
    std::array<double, 3> x_in = {4.3, 2.3, 5.7};
    v1.add_gaussian(x_in);
    double tot_mass = v1.getTot();
    auto tot_com = v1.getWeightedCOM();

    double diff = fabs(1 - v1.getTot());
    double diff2 = 0.0;
    
    for(int i = 0; i < 3; i++){
        diff2 += fabs(x_in[i] - tot_com[i]); 
    }

    if(diff > testing_threshold) {
        std::cout << "Total mass is not being conserved." << std::endl;
        std::cout << "System mass = " << tot_mass << std::endl;
        std::cout << "Expected mass = 1.0" << std::endl;
        return 1;
    }
    if(diff2 > testing_threshold){
        std::cout << "Center of mass is not being conserved." << std::endl;
        for(int i = 0; i < 3; i++){
            std::cout << "Observed COM for dim " << i << " = " << tot_com[i] << std::endl;
            std::cout << "Expected COM for dim " << i << " = " << x_in[i] << std::endl;
        }
        return 1;
    }
    return 0;
}