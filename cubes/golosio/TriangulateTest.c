#include <stdio.h>
#include <stdlib.h>
#include "Triangulate.h"

///////////////////////////////////////////////////////////////////////////
// Iso-Surface Triangulation algorithm:
// test of the version with input/output arrays in memory
// Main function of the test program
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
  float *vol_data, *vertex_data;
  int *triangle_data;
  int nvtx, ntri;
  FILE *vol_fp, *vertex_fp, *triangle_fp;

  if (argc != 12) {
    printf("Usage: triangulate_test inputfile vertexfile trianglefile Nx Ny Nz Sidex Sidey Sidez thresh sign\n");
    exit(0);
  }
  
// read command-line arguments:
  sscanf(argv[4], "%d", &n[0]);
  sscanf(argv[5], "%d", &n[1]);
  sscanf(argv[6], "%d", &n[2]);
  sscanf(argv[7], "%f", &s[0]);
  sscanf(argv[8], "%f", &s[1]);
  sscanf(argv[9], "%f", &s[2]);
  sscanf(argv[10], "%f", &th);
  sscanf(argv[11], "%d", &cmp);

  // allocate memory for volumetric data array
  vol_data = (float*)malloc(sizeof(float)*n[0]*n[1]*n[2]);
  if (!vol_data) {
    printf("Memory allocation error for vol_data in main()");
    exit(EXIT_FAILURE);
  }

  //read volumentric data array from file
  vol_fp = fopen(argv[1], "rb");
  if (!fread (vol_data, sizeof(float), n[0]*n[1]*n[2], vol_fp)) {
    printf("Error reading file in main()\n");
    exit(EXIT_FAILURE);
  }
  fclose(vol_fp);

  // call iso-surface triangulation function
  Triangulate(vol_data, &vertex_data, &triangle_data, n, s, &nvtx, &ntri, th,
	      cmp);

  // write vertex coordinates and gradient vector components output to files
  vertex_fp = fopen(argv[2], "wb");
  fwrite(&nvtx, sizeof(int), 1, vertex_fp);
  fwrite(vertex_data, sizeof(float), nvtx*6, vertex_fp);
  fclose(vertex_fp);

  // write triangle vertex indexes to file
  triangle_fp = fopen(argv[3], "wb");
  fwrite(&ntri, sizeof(int), 1, triangle_fp);
  fwrite(triangle_data, sizeof(int), ntri*3, triangle_fp);
  fclose(triangle_fp);

  //deallocate memory
  free(vol_data);
  free(vertex_data);
  free(triangle_data);

  return 0;
}
