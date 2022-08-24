#include <array>
#include <vector>

template <class T>
using Vec = std::vector<T>;
template <class T>
using Vec2 = std::array<T,2>;
template <class T>
using Vec3 = std::array<T,3>;
struct Tile2d{
  Vec<Vec2<double> > positions;
  Vec2<double> spacing;
  Vec<Vec2<double> > generate_offset_positions(Vec2<int> nx){
    Vec<Vec2<double> > output(positions.size());
    for(int i = 0; i < positions.size(); i++){
      for(int j = 0; j < 2; j++){
        output[i][j] = (positions[i][j] + nx[j])*spacing[j];
      }
    }
    return output;
  }
};

