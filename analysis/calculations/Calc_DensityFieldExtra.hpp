#include "Calc_DensityField.hpp"
#include "../../tools/stlmath.hpp"
#include "../../interface/interface.hpp"
//scheme clusters, omitting grid densities below thresh
//each cluster is averaged, with the sum needing t o

class Calc_DensityFieldExtra : public Calc_DensityField{
  public:
    Calc_DensityFieldExtra(InputPack& input):Calc_DensityField{input}{
      input.params().readString("gro_out", KeyType::Required, groname_);
      input.params().readNumber("starting_density_threshold", KeyType::Required, start_thresh_);
      input.params().readNumber("cluster_density_threshold", KeyType::Required, thresh_);
      input.params().readNumber("span", KeyType::Required, span_);
      input.params().readNumber("mass_cutoff", KeyType::Required, dclust_);      
      return;
    }
    void addNeighbors(const Vec3<int>& coord, std::vector<int>& cluster, std::vector<bool>& isClaimed){
      Vec3<int> coord2 = coord;
      for(int i = -span_; i<=span_; i++){
        for(int j = -span_; j<=span_; j++){
          for(int k = -span_; k<=span_; k++){
            if(i == 0 && j == 0 && k == 0) continue;
            coord2[0] = coord[0]+i; coord2[1] = coord[1]+j; coord2[2] = coord[2]+k;
            if(hasBoxVec_ && (coord2[0] < 0 || coord2[1] < 0 || coord2[2] < 0 || 
                coord2[0] >= npoints_[0] || coord2[1] >= npoints_[1] || coord2[2] >= npoints_[2])) continue;
            int idx = _map31(coord2);
            if(avggridvals_[idx] > thresh_ && !isClaimed[idx]){
              cluster.push_back(idx);
              isClaimed[idx] = 1;
            }
          }        
        }        
      }
      return;
    }
    virtual void finalOutput(){
      Calc_DensityField::finalOutput();
      std::vector<bool> isClaimed(avggridvals_.size(), 0);
      for(int i = 0; i < avggridvals_.size(); i++){
        if(avggridvals_[i] > start_thresh_ && !isClaimed[i]){
          std::vector<int> cluster; 
          if(!isClaimed[i]){
            cluster.push_back(i);
            isClaimed[i] = 1;
          };
          int position = 0;
          while(position < cluster.size()){
            auto coord = _map13(i);
            addNeighbors(coord, cluster, isClaimed);
            position++;
          }
          clusters_.push_back(cluster);
        }
      }
      //now have a list of clusters consisting of indices, need to remove any that don't add up near 1
      for(int i = clusters_.size()-1; i >= 0; i--){
        auto cluster = clusters_[i];
        double sum = 0.0;
        for(auto idx : cluster){
          sum += avggridvals_[idx];
        }
        if(sum > dclust_){
          clusters_.erase(clusters_.begin() + i);
        }
      }
      //now compute centers of mass for complete clusters
      std::vector<Vec3<double> > coms;
      std::vector<Atom> atoms;
      for(auto& cluster : clusters_){
        Vec3<double> com; com.fill(0.0);
        double wtot = 0.0;
        for(auto idx : cluster){
          double weight = avggridvals_[idx];
          auto pos = _map13(idx);
          for(int i = 0; i < 3; i++){
            com[i] += weight*(pos[i]*gridspacing_[i] + minx_[i] + 0.5*gridspacing_[i]);
          }
          wtot += weight;
        }
        com = (1.0/wtot)*com;
        coms.push_back(com);
      }
      Box box_out_;
      box_out_.hasNamedAtoms = 1;
      box_out_.boxvec = box->boxvec;
      int index_iterator = 1;
      for(auto& com : coms){
        Atom a1;
        a1.name = "C"; a1.type = "C"; a1.x = com; a1.v = {0,0,0};
        a1.index = index_iterator; a1.resnr = 1; a1.resname = "MOL";
        box_out_.atoms.push_back(a1);
        index_iterator++; 
      }
      writeGRO(groname_, &box_out_);
      return;
    }
  protected:
    std::string groname_;
    std::vector< std::vector<int> > clusters_;
    //how high a peak has to be to be used as a minima
    double thresh_, dclust_, start_thresh_;
    int span_;
};