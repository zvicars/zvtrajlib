/////////////////////////////////////////
//          OutputData.c               //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
#define OUTPUT_INCREMENT 65536
//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TriangulateGlobal.h"
//---------------------------------------------------------------------------

int Normalize(float *Nrm) { // normalize to 1 the 3D vector Nrm
  int i;
  float n;
  
  n = 0;
  for (i=0; i<3; i++) n += Nrm[i]*Nrm[i];
  n = sqrt(n);
  if (n != 0) {
    for (i=0; i<3; i++) Nrm[i] = Nrm[i]/n;
  }

  return 0;
}

int StoreTriangle(int *IVtx)
//////////////////////////////////////////////////////////////////////
// Store new triangle in file or in memory
////// global variable FileFlag=1 triangle data must be saved in a file
////// global variable FileFlag=0 triangle data must be saved in memory
// RETURN VALUE: zero on success 
// INPUT ARGUMENTS:
// IVtx[3]: indexes of the 3 triangle vertexes 
// OUTPUT:
// If FileFlag=1 the 3 indexes are written to the file pointed by the global
// file pointer TriangleFP as three consecutive integers (int type)
// in binary format.
// If FileFlag=0 the 3 indexes are written to the global array TriangleData
//
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

  if (!FileFlag) { // check if triangle/vertex data must be saved in a file
    if (NTri >= NTriLimit) { // check if it is necessary to reallocate memory
      NTriLimit += OUTPUT_INCREMENT; // increase limit and reallocate memory
      TriangleData = (int*)realloc(TriangleData, sizeof(int)*3*NTriLimit);
      if (!TriangleData) {
	printf("Memory allocation error in StoreTriangle()\n");
	exit(EXIT_FAILURE);
      }
    }
    ipt = 3*NTri;
    for (i=0; i<3; i++) TriangleData[ipt++] = IVtx[i];
    // store the indexes of the 3 triangle vertexes in memory 
  } 
  else fwrite (IVtx, sizeof(int), 3, TriangleFP);
    // write the indexes of the 3 triangle vertexes in the triangle file 
  NTri++;

  return 0;
}

int StoreVertex(float *x, float *nrm)
//////////////////////////////////////////////////////////////////////
// Stores the coordinates and gradient components of a new mesh vertex
// in a file or in memory
////// global variable FileFlag=1 vertex data must be saved in a file
////// global variable FileFlag=0 vertex data must be saved in memory
// RETURN VALUE: zero on success 
// INPUT ARGUMENTS:
// x[3]: vertex x,y,z coordinates
// nrm[3]: gradient x,y,z components
// OUTPUT:
// If FileFlag=1 the vertex coordinates and gradient components are written
// to the file pointed by the global file pointer VertexFP as 6 consecutive
// floating point numbers (float type) in binary format.
// If FileFlag=0 they are written to the global array VertexData
//
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

  if (!FileFlag) { // check if triangle/vertex data must be saved in a file
    if (NVtx >= NVtxLimit) { // check if it is necessary to reallocate memory
      NVtxLimit += OUTPUT_INCREMENT; // increase limit and reallocate memory
      VertexData = (float*)realloc(VertexData, sizeof(float)*6*NVtxLimit);
      if (!VertexData) {
	printf("Memory allocation error in StoreVertex()\n");
	exit(EXIT_FAILURE);
      }
    }
    ipt = 6*NVtx;
    for (i=0; i<3; i++) VertexData[ipt++] = x[i];
    // store the vertex coordinates in memory
    for (i=0; i<3; i++) VertexData[ipt++] = nrm[i];
    // store the 3 gradient vector components in memory 
  } 
  else {
    fwrite (x, sizeof(float), 3, VertexFP);
    // write the vertex coordinates in the vertex file
    fwrite (nrm, sizeof(float), 3, VertexFP);
    // write the 3 gradient vector components in the vertex file 
  }

  NVtx++;

  return NVtx-1;
}

int OutputDataInit() // output data initialization
{

  NVtxLimit = 0;
  NTriLimit = 0;

  return 0;
}

int OutputDataFree() // free the memory allocated for output data
{

  free(VertexData);
  free(TriangleData);

  return 0;
}




