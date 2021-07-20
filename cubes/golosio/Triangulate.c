/////////////////////////////////////////
//          Triangulate.c              //
//     Author : Bruno Golosio          //
/////////////////////////////////////////

//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "TriangulateGlobal.h"
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
int Triangulate(float *vol_data, float **vertex_data,  int **triangle_data,
		int *n, float *s, int *nvtx, int *ntri, float mu, int cmp)
{
  int i, iz;
  
  FileFlag = 0; // output array are stored in memory rather than in files
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
  }

  TriangulateFree(); // deallocate memory

  // set output pointers 
  *vertex_data = VertexData;
  *triangle_data = TriangleData;
  *nvtx = NVtx;
  *ntri = NTri;

  return 0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// TriangulateFile: Iso-surface triangulation algorithm
// (version with input volumetric data from file, output stored in files)
// INPUT:
// InputFile : 3D volumetric data file containing Nx*Ny*Nz elements
////// binary (raw) format:
////// Nx*Ny*Nz floats, Nx running faster
// n : array dimensions Nx=n[0], Ny=n[1], Nz=n[2]
// s : voxel sizes sx=s[0], sy=s[1], sz=s[2]
// mu : Iso-Surface Threshold Value (iso-value)
// cmp : Sign defining volume inside iso-surface
//       (< , <=, >=, >) : LT_THRESH, LE_THRESH, GE_THRESH, GT_THRESH
// OUTPUT:
// ntri: number of triangles
// nvtx: number of vertexes
// VertexFile : file with vertex coordinates and gradient vector components
////// binary format:
////// Nv (int type): number of vertexes in the mesh
////// x1, y1, z1 (float type): coordinates of vertexes n. 1
////// vx1, vy1, vz1 (float type): gradient components of vertexes n. 1
////// ............
////// xNv, yNv, zNv (float type): coordinates of vertexes n. Nv
////// vxNv, vyNv, vzNv (float type): gradient components of vertexes n. Nv
// TriangleFile : file with triangle vertex indexes
////// binary format:
////// Nt (int type): number of triangles in the mesh
////// i1, j1, k1 (float type): indexes of the mesh vertexes corresponding 
////// to the 3 vertexes of triangle n. 1
////// ............
////// iNt, jNt, kNt (float type): indexes of the mesh vertexes corresponding 
////// to the 3 vertexes of triangle n. Nt
//
// DESCRIPTION:
// Iso-surface triangulation algorithm (version with input volumetric data
// from file, output stored in files). It makes a loop on all slices
// of the input volumetric data array vol_data. The slice index is iz.
// At each step of the loop, it keeps 4 consecutive slices in memory
// by shifting the last 3 slices and loading from the file
// InputFile a new one through the
// function  LoadSlices(iz). Then it evaluates volumetric data gradients
// at the logical cube edges through the function MkGrad(iz),
// it evaluates intersection of isosurface with logical cube edges
// through the function MkEdge(iz), and it perform a standard triangulation 
// inside each logical cube through the function PolygonizeCube(iz)
// For each vertex of the resulting mesh, the 3 vertex coordinates and the
// 3 gradient vector components are written to the file VertexFile
// as 6 consecutive floating point numbers (float type) in binary format.
// For each triangle of the resulting mesh, the 3 indexes associated to the
// triangle vertexes are written to the file TriangleFile
// as 3 consecutive integers (int type) in binary format
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int TriangulateFile(char *InputFile, char *VertexFile,  char *TriangleFile,
		int *n, float *s, float mu, int cmp)
{
  int i, iz;
  
  FileFlag = 1; // output array are stored in files rather than in memory
  // initializations:
  NVtx = 0;
  NTri = 0;
  IsoValue = mu;
  CmpSign = cmp;
  for (i=0; i<3; i++) {
    NVxl[i] = n[i];
    Side[i] = s[i];
  }
  
  // open and initializes input/output files
  VoxelFP = fopen(InputFile, "rb");
  VertexFP = fopen(VertexFile, "wb");
  fwrite (&NVtx, sizeof(int), 1, VertexFP);
     // leaves a space for the number of vertexes
  TriangleFP = fopen(TriangleFile, "wb");
  fwrite (&NTri, sizeof(int), 1, TriangleFP);
     // leaves a space for the number of triangles
  TriangulateInit();
  
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
  fseek(VertexFP, 0, SEEK_SET);
  fwrite (&NVtx, sizeof(int), 1, VertexFP); // write number of vertexes
  fseek(TriangleFP, 0, SEEK_SET);
  fwrite (&NTri, sizeof(int), 1, TriangleFP); // write number of triangles

  // close files
  fclose(VoxelFP);
  fclose(VertexFP);
  fclose(TriangleFP);
  TriangulateFree(); // deallocate memory

  return 0;
}

//---------------------------------------------------------------------------
int LoadSlices(int iz)    // load 4 slices in memory
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
    if (!FileFlag) {
      memcpy(&Vertex[3][0][0], &VolData[NVxl[0]*NVxl[1]*(iz+2)],
	     sizeof(float)*NVxl[0]*NVxl[1]);
    }
    else {
      if (!fread (&Vertex[3][0][0], sizeof(float), NVxl[0]*NVxl[1], VoxelFP)) {
	printf("Error reading file in LoadSlices().\n");
	// printf("%d\t%d\t%d\n", ix, iy, iz);
	exit(EXIT_FAILURE);
      }
    }
  }

  return 0;
}

int TriangulateInit() // Initializes arrays used by the program
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

int TriangulateFree() // deallocate memory of arrays used by the program
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
