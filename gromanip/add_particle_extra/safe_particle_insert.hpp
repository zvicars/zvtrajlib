//header file for an algorithm to inexpensively create random particle packings within a volume
//general idea is that there's a 3D matrix of "safe" gridcells to place an atom
//particles will be placed randomly until it's difficult to find random spots to place them
//then an explicit list of safe cells will be generated for placing the remainging atoms
#pragma once
#include "../volume.hpp"
#include "../../tools/Matrix.hpp"
#include <random>
class safeInsert{
    public:
    safeInsert(double gridspacing, Volume* volume_ptr){
        //xmin ymin zmin xmax ymax zmax
        auto bounding_box = volume_ptr->getBoundingBox();
        gridspacing_=gridspacing;
        for(int i=0; i<3; i++){
            box_min_[i] = bounding_box[i];
            box_max_[i] = bounding_box[i+3];
        }

        Vec3<std::size_t> ncells;
        for(int i=0; i<3; i++){
            ncells[i] = (box_max_[i] - box_min_[i])/gridspacing_;
        }
        field_.initialize(ncells);
        field_.fill(0);
        for(int i = 0; i < field_.size_1d(); i++){
            auto index = field_.map1N(i);
            Vec3<double> real_pos;
            for(int j = 0; j < 3; j++){
                real_pos[j] = index[j]*gridspacing_ +  box_min_[j];
            }
            field_.at_1d(i) = (int)!volume_ptr->isInside(real_pos);
        }
        return;
    }
    void fillAroundIndex(Vec3<std::size_t> index, double radius){
        double modified_radius = radius + std::sqrt(3)*gridspacing_;
        int bound = std::ceil(modified_radius / gridspacing_);
        Vec3<double> real_pos;
        for(int l = 0; l<3; l++){
          real_pos[l] = index[l]*gridspacing_ + box_min_[l];
        }
        auto box_ncells = field_.size();
        for(int i = (int)index[0] - bound; i <= (int)index[0] + bound; i++){
          for(int j = (int)index[1] - bound; j <= (int)index[1] + bound; j++){
            for(int k = (int)index[2] - bound; k <= (int)index[2] + bound; k++){
              if(i < 0 || i >= box_ncells[0] || j < 0 || j >= box_ncells[1] || k < 0 || k >= box_ncells[2]){
                continue;
              }
              Vec3<std::size_t> index2 = {(std::size_t)i, (std::size_t)j, (std::size_t)k};
              Vec3<double> real_pos_test;
              for(int l = 0; l<3; l++){
                real_pos_test[l] = index2[l]*gridspacing_ + box_min_[l];
              }
              real_pos_test = {i*gridspacing_ + box_min_[0], j*gridspacing_ + box_min_[1], k*gridspacing_ + box_min_[2]};
              double dist = norm2(real_pos_test - real_pos);
              if(dist < modified_radius) field_.at(index2) = 1;
            }
          }
        }
        return;
    }
    
    Vec3<std::size_t> getSafeIndex(){
        //rng_policy = 0 implies it's going to be randomly searching for sites, rng_policy = 1 involves
        //explicit site search 
        if(rng_policy == 0){
            std::uniform_int_distribution<int> distribution(0, field_.size_1d());
            for(int i = 0; i < 10; i++){
                std::size_t test_index = distribution(generator);
                if(field_.at_1d(test_index) == 0) return field_.map1N(test_index);
            }
            std::cout << "changing rng policy" << std::endl;
            rng_policy = 1;
        }
        return getSafeIndexVector();
    }
    Vec3<std::size_t> getSafeIndexVector(){
        safe_indices_.reserve(field_.size_1d());
        //if the safe-indiced vector hasn't been built, build it
        if(!sivInitialized){
            for(int i = 0; i < field_.size_1d(); i++){
                if(field_.at_1d(i) == 0){
                    safe_indices_.push_back(i);
                }
            }
            sivInitialized=1;
        }
        //otherwise, make sure the entries are still safe
        else{
            for(int i=0; i < safe_indices_.size(); i++){
                if(field_.at_1d(safe_indices_[i]) == 1){
                    safe_indices_.erase(safe_indices_.begin() + i);
                    i--;
                }
            }
        }
        std::size_t size = safe_indices_.size();
        if(size == 0){
            fail_ = 1;
            return {0,0,0};
        }
        std::uniform_int_distribution<int> distribution(0, size);
        std::size_t index = safe_indices_[distribution(generator)];
        return field_.map1N(index);
    }
    bool getSafePosition(Vec3<double>& pos, double radius){
        auto index = getSafeIndex();
        /*std::ofstream ofile("test1.txt");
        for(int i = 0; i < field_.size_1d(); i++){
            auto idx3d = field_.map1N(i);
            if(field_.at_1d(i) == 0) ofile << idx3d[0] << " " << idx3d[1] << " " << idx3d[2] << "\n";
        }
        ofile << std::endl;
        ofile.close();
        std::cin.get();*/
        fillAroundIndex(index, radius);
        /*ofile.open("test2.txt");
        for(int i = 0; i < field_.size_1d(); i++){
            auto idx3d = field_.map1N(i);
            if(field_.at_1d(i) == 0) ofile << idx3d[0] << " " << idx3d[1] << " " << idx3d[2] << "\n";
        }
        ofile << std::endl;
        ofile.close();
        std::cin.get(); */
        for(int i = 0; i < 3; i++){
            pos[i] = index[i]*gridspacing_ + box_min_[i];
        }
        return fail_;
    }
    private:
    Matrix<int,3> field_;
    Vec3<double> box_min_;
    Vec3<double> box_max_;
    double gridspacing_;
    std::vector<std::size_t> safe_indices_;
    //rng stuff
    bool fail_ = 0; //have I run out of safe sites?
    std::default_random_engine generator;
    int rng_policy=0;
    bool sivInitialized=0;
};