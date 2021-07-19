/////////////////////////////////////////
//          Triangulate.h              //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
#ifndef TriangulateH
#define TriangulateH

#include <stdio.h>

// Triangulate: Iso-surface triangulation algorithm
// (version with input/output array in memory)
int Triangulate(float *vol_data, float **vertex_data,  int **triangle_data,
		int *n, float *s, int *nvtx, int *ntri, float thresh, int cmp);

// TriangulateFile: Iso-surface triangulation algorithm
// (version with input volumetric data from file, output stored in files)
int TriangulateFile(char *InputFile, char *VertexFile,  char *TriangleFile,
		int *n, float *s, float thresh, int cmp);

#endif

