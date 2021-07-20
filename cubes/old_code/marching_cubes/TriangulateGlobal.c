#include <stdio.h>

FILE *VertexFP;  // file of vertexes
FILE *TriangleFP;  // file of triangles
FILE *VoxelFP; // file of voxels
int NVtx;   // n. of vertexes
int NTri;  // n. of triangles
int NVtxLimit;   // max n. of vertexes (dinamically updated)
int NTriLimit;  // max n. of triangles (dinamically updated)
int NVxl[3]; // voxel array size
float Side[3]; // voxel sides (x, y, z)
float IsoValue; // Iso-Surface Threshold Value
int CmpSign; // Sign defining volume inside iso-surface
             // LT_THRESH, LE_THRESH, GE_THRESH, GT_THRESH 
float **Vertex[4];  // Mu Values in 4 horizontal planes of cube vertexes
float **Grad[3][2]; // gradient in cube vertexes [component][z plane][y][x]
int **EdgeIpt[3][2]; // indexes of intersection points of the IsoSurface
                     // with the cube edges [edge axis][z plane][y][x]
float *VolData; // array of volumetric data
float *VertexData; // array of vertexes
int *TriangleData; // array of triangles
int FileFlag; // input/output flag (0: from/to arrays, 1: from/to files)
