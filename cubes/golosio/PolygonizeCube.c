/////////////////////////////////////////
//          PolygonizeCube.c           //
//     Author : Bruno Golosio          //
/////////////////////////////////////////

//---------------------------------------------------------------------------
#include "TriangulateGlobal.h"
#include "TriangulateTab.h"
//---------------------------------------------------------------------------

int dx[8] = {1, 1, 1, 1, 0, 0, 0, 0}; // relative position of cube vertexes
int dy[8] = {0, 1, 1, 0, 0, 1, 1, 0};
int dz[8] = {0, 0, 1, 1, 0, 0, 1, 1};
int Edgx[12] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}; //
int Edgy[12] = {0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1};
int Edgz[12] = {0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1};
int Edga[12] = {1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 0, 0};

int PolygonizeCube(int iz)
//////////////////////////////////////////////////////////////////////
// Performs a standard triangulation inside each logical cube
// RETURN VALUE: zero on success 
// INPUT ARGUMENTS:
// - iz: index of slice from volumetric data
// INPUT FROM GLOBAL VARIABLES:
// **Vertex[4]: Mu Values of cube vertexes in 4 slices
// **EdgeIpt[3][2]: indexes of intersection points of the IsoSurface
//                  with the cube edges [edge axis][z plane][y][x]
// alntri[256]: (from lookup table) alntri[ib] is the number of triangles
//              inside logical cube for the standard triangulation num. ib 
// alipt[256][5][3]: (from lookup table) alipt[ib][itr][iv]
//      for the standard triangulation num. ib, for the triangle num. itr
//      of that triangulation, and for vertex num. iv of that triangle,
//      alipt[ib][itr][iv] is the index of the cube edge on which the vertex
//      is positioned
// OUTPUT TO FILE (OR MEMORY):
// For each triangle of the mesh, the 3 indexes associated to the triangle
// vertexes are stored in a file (or in memory) through a call to the function
// "StoreTriangle"
//
// DESCRIPTION:
// This function is called by the function "Triangulate" (or "TriangulateFile")
// as part of the marching cube algorithm.
// It makes a loop on the logical cubes of a volumetric data slice. For
// each logical cube it checks if the 8 cube vertexes are above or below the
// isovalue, thus obtaining 8 boolean values. Those values are converted to
// a binary number, than the standard triangulation corresponding to that
// binary number is retrieved from the lookup table.
// This standard triangulation says how many triangles are present inside
// the logical cube and on which edges of the logical cube the triangle
// vertexes are positioned. In this way the algorithm associates to each
// triangle vertex the index of the corresponding mesh vertex.
// The 3 indexes associated to the triangle vertexes are stored in a file
// (or in memory) through a call to the function "StoreTriangle"
//////////////////////////////////////////////////////////////////////
{

  int InVert[8]; // boolean array of cube vertexes
                 // (true, false) = (inside, outside)
  int ix, iy;
  int i;
  float Mu;
  int VertexByte;   // byte from InVert[]
  int itr, jpt; // indexes for triangles, vertexes
  int VtxIdx[3];
  int iedg;
  
  for (iy=0; iy<(NVxl[1]-1); iy++) {     // loop on cubes
  for (ix=0; ix<(NVxl[0]-1); ix++) {
    for (i=0; i<8; i++) {          // loop on cube vertexes index
      Mu = Vertex[dz[i]][iy+dy[i]][ix+dx[i]]; // Mu Value of i-th vertex
      if ((CmpSign==LT_THRESH && Mu<IsoValue) ||
	  (CmpSign==LE_THRESH && Mu<=IsoValue) ||
	  (CmpSign==GE_THRESH && Mu>=IsoValue) ||
	  (CmpSign==GT_THRESH && Mu>IsoValue)) 
	InVert[i] = 1;             // compares with IsoValue
      else InVert[i] = 0;
    }
    VertexByte = BoolToByte(InVert); // converts InVert[] into a byte
    // VertexByte is an index that points to the standard triangulation
    // to be performed within the cube

    for (itr=0; itr<alntri[VertexByte]; itr++) { // loop on triangles
      for (jpt=0; jpt<3; jpt++) {  // loop on triangle vertexes
	iedg = alipt[VertexByte][itr][jpt];
	// index of the cube edge on which the vertex is positioned
	VtxIdx[jpt] = EdgeIpt[Edga[iedg]][Edgz[iedg]]
	  [iy+Edgy[iedg]][ix+Edgx[iedg]]; // index of the point
      }
      StoreTriangle(VtxIdx); // store the triangle
    }
  }
  }

  return 0;
}

// converts a boolean array InVert[8] to a byte VertexByte
int BoolToByte(int InVert[]) {
  int VertexByte;
  int i;
  
  VertexByte = 0;
  for (i=7; i>=0; i--) {
    VertexByte <<= 1;
    if(InVert[i]) VertexByte ++;
  }
  return VertexByte;
}
