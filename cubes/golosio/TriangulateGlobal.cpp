#include "TriangulateGlobal.hpp"
#include <cmath>
#include <cstring>
int TriangulateGlobal::Triangulate(float *vol_data, std::vector<float>& vertex_data,  std::vector<int>& triangle_data, int *n, float *s, int *nvtx, int *ntri, float mu, int cmp)
{
  int i, iz;
  // initializations:
  NVtx = 0;
  NTri = 0;
  IsoValue = mu;
  CmpSign = cmp;
  for (i=0; i<3; i++) {
    NVxl[i] = n[i];
    Side[i] = s[i];
  }
  
  VolData = vol_data;
  TriangulateInit();
  OutputDataInit();
  
  LoadSlices(-2);                   // load 4 slices in memory
  LoadSlices(-1);                   // load 4 slices in memory
  MkGrad(-1);                       // evaluate gradient
  for (iz=0; iz<NVxl[2]; iz++) {
    LoadSlices(iz);                   // load 4 slices in memory
    MkGrad(iz);                       // evaluate gradient
    MkEdge(iz);  // evaluate intersection of isosurface with cube edges
    if (iz > 0) PolygonizeCube(iz); // perform standard triangulations
    //printf("%d / %d\n", iz+1, NVxl[2]);   // update progress bar
  }

  TriangulateFree(); // deallocate memory

  // set output pointers 
  vertex_data = VertexData;
  triangle_data = TriangleData;
  *nvtx = NVtx;
  *ntri = NTri;

  return 0;
}

int TriangulateGlobal::LoadSlices(int iz)    // load 4 slices in memory
{
  int ix, iy;
  for (iy=0; iy<NVxl[1]; iy++) {
  for (ix=0; ix<NVxl[0]; ix++) {
    // The first three slices are updated through a scroll
    Vertex[0][iy][ix] = Vertex[1][iy][ix];
    Vertex[1][iy][ix] = Vertex[2][iy][ix];
    Vertex[2][iy][ix] = Vertex[3][iy][ix];
    // the 4-th slice is loaded from the voxel file or array
  }
  }
  if (iz<(NVxl[2]-2)) {
    memcpy(&Vertex[3][0][0], &VolData[NVxl[0]*NVxl[1]*(iz+2)],
    sizeof(float)*NVxl[0]*NVxl[1]);
  }
  return 0;
}

int TriangulateGlobal::TriangulateInit() // Initializes arrays used by the program
{
  int i, j, ix, iy;
  for (i=0; i<4; i++) Vertex[i] = float_matrix(NVxl[1], NVxl[0]);
  for (j=0; j<2; j++) {
    for (i=0; i<3; i++) {
      EdgeIpt[i][j] = int_matrix(NVxl[1], NVxl[0]);
      Grad[i][j] = float_matrix(NVxl[1], NVxl[0]);
    }
  }
  for (iy=0; iy<NVxl[1]; iy++) {
    for (ix=0; ix<NVxl[0]; ix++) {
      Vertex[1][iy][ix] = 0;
      Vertex[2][iy][ix] = 0;
      Vertex[3][iy][ix] = 0;
      EdgeIpt[0][1][iy][ix] = 0;
      EdgeIpt[1][1][iy][ix] = 0;
      Grad[0][1][iy][ix] = 0;
      Grad[1][1][iy][ix] = 0;
      Grad[2][1][iy][ix] = 0;
    }
  }

 return 0;
}


int TriangulateGlobal::TriangulateFree() // deallocate memory of arrays used by the program
{
 int i, j;

 for (i=0; i<4; i++) free_float_matrix(Vertex[i]);
 for (j=0; j<2; j++) {
   // free_int_matrix(CubeFlag[j]);
   for (i=0; i<3; i++) {
     free_int_matrix(EdgeIpt[i][j]);
     free_float_matrix(Grad[i][j]);
   }
 }
 return 0;
}

int TriangulateGlobal::StoreTriangle(int *IVtx)
//////////////////////////////////////////////////////////////////////
// DESCRIPTION:
// This function is called by the function "PolygonizeCube"
// as part of the marching cube algorithm.
// It stores a new triangle in a file (or in memory, depending on the flag
//  FileFlag) by writing 3 indexes associated to the triangle vertexes
// in the file (or in memory) as 3 consecutive integers (int type) in
// binary format
//////////////////////////////////////////////////////////////////////
{
  int i, ipt;
  if (NTri >= NTriLimit) { // check if it is necessary to reallocate memory
    NTriLimit += OUTPUT_INCREMENT; // increase limit and reallocate memory
    TriangleData.resize(3*NTriLimit);
  }
  ipt = 3*NTri;
  for (i=0; i<3; i++) TriangleData[ipt++] = IVtx[i];
  // store the indexes of the 3 triangle vertexes in memory 
  // write the indexes of the 3 triangle vertexes in the triangle file 
  NTri++;

  return 0;
}

int TriangulateGlobal::StoreVertex(float *x, float *nrm)
//////////////////////////////////////////////////////////////////////
// DESCRIPTION:
// This function is called by the function "MkEdge"
// as part of the marching cube algorithm.
// It stores the coordinates and gradient components of a new mesh vertex
// in a file (or in memory, depending on the flag
//  FileFlag) by writing the 3 coordinates and the 3 gradient components
// in the file (or in memory) as 6 consecutive floating point numbers
// (float type) in binary format
//////////////////////////////////////////////////////////////////////
{
  int i, ipt;
  if (NVtx >= NVtxLimit) { // check if it is necessary to reallocate memory
    NVtxLimit += OUTPUT_INCREMENT; // increase limit and reallocate memory
    //VertexData = (float*)realloc(VertexData, sizeof(float)*6*NVtxLimit);
    VertexData.resize(6*NVtxLimit);
  }
  ipt = 6*NVtx;
  for (i=0; i<3; i++) VertexData[ipt++] = x[i];
  // store the vertex coordinates in memory
  for (i=0; i<3; i++) VertexData[ipt++] = nrm[i];
  // store the 3 gradient vector components in memory 
  NVtx++;
  return NVtx-1;
}

int TriangulateGlobal::OutputDataInit() // output data initialization
{

  NVtxLimit = 0;
  NTriLimit = 0;

  return 0;
}

int TriangulateGlobal::OutputDataFree() // free the memory allocated for output data
{
  return 0;
}