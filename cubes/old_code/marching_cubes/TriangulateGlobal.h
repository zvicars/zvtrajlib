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
//---------------------------------------------------------------------------
// Functions
//---------------------------------------------------------------------------
int LoadSlices(int iz);         // load 4 slices in memory
int TriangulateInit(); // initialize arrays
int TriangulateFree(); // free arrays
int MkEdge(int iz);  // evaluate intersection points of the IsoSurface
                      // with the logical cubes edges and volumetric data
                      // gradient components in such points
int MkGrad(int iz);             // evaluate gradients in logical cube vertexes
int PolygonizeCube(int iz); // Performs a standard triangulation
                            // in each logical cube
int StoreTriangle(int VtxIdx[]); // stores triangles in a file (or in memory)
int BoolToByte(int InVert[]);   // converts an array of 8 booleans into a byte
int StoreVertex(float *x, float *nrm); // converts an array of 8 booleans
                                       // into a byte
int OutputDataInit(); // output data initialization
int OutputDataFree(); // output data deallocate memory

//---------------------------------------------------------------------------
// Variables
//---------------------------------------------------------------------------
extern FILE *VertexFP;    // file of vertexes
extern FILE *TriangleFP;  // file of triangles
extern FILE *VoxelFP; // file of voxels
extern int NVtx;   // n. of vertexes
extern int NTri;  // n. of triangles
extern int NVtxLimit;   // max n. of vertexes (dinamically updated)
extern int NTriLimit;  // max n. of triangles (dinamically updated)
extern int NVxl[3]; // voxel array size
extern float Side[3];        // voxel sides (x, y, z)
extern float IsoValue; // Iso-Surface Threshold Value
extern int CmpSign; // Sign defining volume inside iso-surface
                    // LT_THRESH, LE_THRESH, GE_THRESH, GT_THRESH 
extern float **Vertex[4];  // Mu Values of logical cube vertexes in 4 slices
extern float **Grad[3][2]; // gradient in logical cube vertexes
                           // [component][z plane][y][x]
extern int **EdgeIpt[3][2]; // indexes of intersection points of the IsoSurface
                            // with the logical cube edges [edge axis][z plane][y][x]
extern float *VolData; // array of volumetric data
extern float *VertexData; // array of vertexes
extern int *TriangleData; // array of triangles
extern int FileFlag; // input/output flag (0: from/to arrays, 1: from/to files)

#endif
