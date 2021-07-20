/////////////////////////////////////////
//          TriangulateGlobal.h        //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
#ifndef TriangulateGlobalH
#define TriangulateGlobalH
#define LT_THRESH 0
#define LE_THRESH 1
#define GE_THRESH 2
#define GT_THRESH 3
#include <stdio.h>
#include "matrix.hpp"
#include <vector>
class TriangulateGlobal{
public:
  int Triangulate(float *vol_data, std::vector<float>& vertex_data,  std::vector<int>& triangle_data, int *n, float *s, int *nvtx, int *ntri, float mu, int cmp);
private:
  int LoadSlices(int iz);         // load 4 slices in memory
  int TriangulateInit(); // initialize arrays
  int TriangulateFree(); // free arrays
  int MkEdge(int iz);  // evaluate intersection points of the IsoSurface
  int MkGrad(int iz);             // evaluate gradients in logical cube vertexes
  int PolygonizeCube(int iz); // Performs a standard triangulation
  int StoreTriangle(int VtxIdx[]); // stores triangles in a file (or in memory)
  int BoolToByte(int InVert[]);   // converts an array of 8 booleans into a byte
  int StoreVertex(float *x, float *nrm); // converts an array of 8 booleans into byte
  int OutputDataInit(); // output data initialization
  int OutputDataFree(); // output data deallocate memory

  //---------------------------------------------------------------------------
  // Variables
  //---------------------------------------------------------------------------
  int NVtx;   // n. of vertexes
  int NTri;  // n. of triangles
  int NVtxLimit;   // max n. of vertexes (dinamically updated)
  int NTriLimit;  // max n. of triangles (dinamically updated)
  int NVxl[3]; // voxel array size
  float Side[3];        // voxel sides (x, y, z)
  float IsoValue; // Iso-Surface Threshold Value
  int CmpSign; // Sign defining volume inside iso-surface
  float **Vertex[4];  // Mu Values of logical cube vertexes in 4 slices
  float **Grad[3][2]; // gradient in logical cube vertexes
  int **EdgeIpt[3][2]; // indexes of intersection points of the IsoSurface
  float *VolData; // array of volumetric data
  std::vector<float> VertexData; // array of vertexes
  std::vector<int> TriangleData; // array of triangles
  int OUTPUT_INCREMENT = 65536;
};
#endif
