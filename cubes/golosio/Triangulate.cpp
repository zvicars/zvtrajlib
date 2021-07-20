/////////////////////////////////////////
//          Triangulate.c              //
//     Author : Bruno Golosio          //
/////////////////////////////////////////

//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TriangulateGlobal.hpp"
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Triangulate: Iso-surface triangulation algorithm
// (version with input/output array in memory)
// INPUT:
// vol_data[Nx*Ny*Nz] : 3D volumetric data array
// n : array dimensions Nx=n[0], Ny=n[1], Nz=n[2]
// s : voxel sizes sx=s[0], sy=s[1], sz=s[2]
// mu : Iso-Surface Threshold Value (iso-value)
// cmp : Sign defining volume inside iso-surface
//       (< , <=, >=, >) : LT_THRESH, LE_THRESH, GE_THRESH, GT_THRESH
// OUTPUT:
// ntri: number of triangles
// nvtx: number of vertexes
// vertex_data : array of vertex coordinates and gradient vector components
// triangle_data : array of triangle vertex indexes
//
// DESCRIPTION:
// Iso-surface triangulation algorithm. It makes a loop on all slices
// of the input volumetric data array vol_data. The slice index is iz.
// At each step of the loop, it keeps 4 consecutive slices in memory
// by shifting the last 3 slices and loading a new one through the
// function  LoadSlices(iz). Then it evaluates volumetric data gradients
// at the logical cube edges through the function MkGrad(iz),
// it evaluates intersection of isosurface with logical cube edges
// through the function MkEdge(iz), and it perform a standard triangulation 
// inside each logical cube through the function PolygonizeCube(iz)
// For each vertex of the resulting mesh, the 3 vertex coordinates and the
// 3 gradient vector components are written to the array vertex_data
// as 6 consecutive floating point numbers (float type).
// For each triangle of the resulting mesh, the 3 indexes associated to the
// triangle vertexes are written to the array triangle_data
// as 3 consecutive integers (int type)
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int TriangulateGlobal::Triangulate(float *vol_data, float **vertex_data,  int **triangle_data,
		int *n, float *s, int *nvtx, int *ntri, float mu, int cmp)
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
    printf("%d / %d\n", iz+1, NVxl[2]);   // update progress bar
  }

  TriangulateFree(); // deallocate memory

  // set output pointers 
  *vertex_data = VertexData;
  *triangle_data = TriangleData;
  *nvtx = NVtx;
  *ntri = NTri;

  return 0;
}