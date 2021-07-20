#include <stdio.h>
#include <stdlib.h>
#include "Triangulate.h"

///////////////////////////////////////////////////////////////////////////
// Iso-Surface Triangulation algorithm
// Main function of the program
// read command-line arguments:
// inputfile : 3D volumetric data file containing Nx*Ny*Nz elements
// vertexfile : file with vertex coordinates and gradient vector components
// trianglefile : file with triangle vertex indexes
// Nx, Ny, Nz : volumetric data array dimensions
// Sidex, Sidey, Sidez : voxel sizes
// thresh : Iso-Surface Threshold Value (iso-value)
// sign : Sign defining volume inside iso-surface
//       (< , <=, >=, >) : LT_THRESH, LE_THRESH, GE_THRESH, GT_THRESH
///////////////////////////////////////////////////////////////////////////
// DESCRIPTION:
// Main function for Iso-surface triangulation using an improved marching cube
// algorithm. It loads volumetric data from the file inputfile, builds a
// triangulated model of the isosurface corresponding to the isovalue thresh,
// and write the resulting mesh in two files, one (vertexfile) containing
// the vertex coordinates and gradient vector components, the other
// (trianglefile) containing, for each triangle of the mesh, the 3
// indexes of the mesh vertexes corresponding to the triangle vertexes 
// INPUT FILE FORMAT:
// inputfile : 3D volumetric data file containing Nx*Ny*Nz elements
////// binary (raw) format:
////// Nx*Ny*Nz floats, Nx running faster
// OUTPUT FILE FORMAT:
// vertexfile : file with vertex coordinates and gradient vector components
////// binary format:
////// Nv (int type): number of vertexes in the mesh
////// x1, y1, z1 (float type): coordinates of vertexes n. 1
////// vx1, vy1, vz1 (float type): gradient components of vertexes n. 1
////// ............
////// xNv, yNv, zNv (float type): coordinates of vertexes n. Nv
////// vxNv, vyNv, vzNv (float type): gradient components of vertexes n. Nv
// trianglefile : file with triangle vertex indexes
////// binary format:
////// Nt (int type): number of triangles in the mesh
////// i1, j1, k1 (float type): indexes of the mesh vertexes corresponding 
////// to the 3 vertexes of triangle n. 1
////// ............
////// iNt, jNt, kNt (float type): indexes of the mesh vertexes corresponding 
////// to the 3 vertexes of triangle n. Nt
//////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  int n[3];
  float s[3];
  float th;
  int cmp;

  if (argc != 12) {
    printf("Usage: triangulate inputfile vertexfile trianglefile Nx Ny Nz Sidex Sidey Sidez thresh sign\n");
    exit(0);
  }

  // read commmand line arguments  
  sscanf(argv[4], "%d", &n[0]);
  sscanf(argv[5], "%d", &n[1]);
  sscanf(argv[6], "%d", &n[2]);
  sscanf(argv[7], "%f", &s[0]);
  sscanf(argv[8], "%f", &s[1]);
  sscanf(argv[9], "%f", &s[2]);
  sscanf(argv[10], "%f", &th);
  sscanf(argv[11], "%d", &cmp);

  TriangulateFile(argv[1], argv[2], argv[3], n, s, th, cmp);

  return 0;
}
