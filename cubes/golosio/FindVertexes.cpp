/////////////////////////////////////////
//          FindVertexes.c             //
//     Author : Bruno Golosio          //
/////////////////////////////////////////

//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include "TriangulateGlobal.hpp"
//---------------------------------------------------------------------------
#define EPS 1e-2

int TriangulateGlobal::MkEdge(int iz)
//////////////////////////////////////////////////////////////////////
// Evaluates intersection points of the IsoSurface with the logical cube edges
// and volumetric data gradiet at the intersection positions.
// RETURN VALUE: zero on success 
// INPUT ARGUMENTS:
// - iz: index of slice from volumetric data
// INPUT FROM GLOBAL VARIABLES:
// **Vertex[4]: Mu Values of logical cube vertexes in 4 slices
// **Grad[3][2]: gradient in logical cube vertexes [component][z plane][y][x]
// OUTPUT TO GLOBAL VARIABLES AND TO FILE (OR MEMORY):
// **EdgeIpt[3][2]: indexes of intersection points of the IsoSurface
//                  with the cube edges [edge axis][z plane][y][x]
// Intersection point coordinates and gradient at the intersection positions
// are stored in a file (or in memory) through a call to the function
// "StoreVertex"
//
// DESCRIPTION:
// This function is called by the function "Triangulate" (or "TriangulateFile")
// as part of the marching cube algorithm.
// It makes a loop on the logical cubes of a volumetric data slice. For
// each vertex of the cube it checks whether it is intersected by the
// isosurface, and in such case it evaluates the intersection position
// and the volumetric data gradient at the intersection position 
// by linear interpolation as described in the article.
// The index of the intersection is stored in the global array EdgeIpt
// The intersection point coordinates and gradient at the intersection
// positions are stored in a file (or in memory) through a call to the
// function "StoreVertex"
//////////////////////////////////////////////////////////////////////
{
  int i, j, dz;
  int ix[3], ix1[3];
  float t, Mu0, Mu1;
  float x[3], nrm[3];
  
  ix[2]=iz;
  for (ix[1]=0; ix[1]<NVxl[1]; ix[1]++) {    // loop on cube vertexes
  for (ix[0]=0; ix[0]<NVxl[0]; ix[0]++) {
    // cube vertex coordinates
      
    Mu0 = Vertex[1][ix[1]][ix[0]];           // Mu of first vertex
    for (j=0; j<3; j++) {   // loop on edge axis (0, 1, 2) = (x, y, z
      dz = (j==2) ?1 : 0;
      EdgeIpt[j][0][ix[1]][ix[0]] = EdgeIpt[j][1][ix[1]][ix[0]];
      // update cube edges plane
      EdgeIpt[j][1][ix[1]][ix[0]] = 0;
      // remains 0 if the IsoSurface does not intersect the edge
      for (i=0; i<3; i++) ix1[i] = ix[i];
      ix1[j]++;
      if (ix1[j] < NVxl[j]) {
	Mu1 = Vertex[1+dz][ix1[1]][ix1[0]]; // Mu value of second vertex
	if (((CmpSign==GT_THRESH || CmpSign==LE_THRESH) && 
	     ((Mu0>IsoValue) != (Mu1>IsoValue))) ||
	    ((CmpSign==GE_THRESH || CmpSign==LT_THRESH) && 
	     ((Mu0>=IsoValue) != (Mu1>=IsoValue)))) {
	  t = (Mu0 - IsoValue)/(Mu0 - Mu1);
	  if (t < EPS) t = EPS;
	  if (t > 1-EPS) t = 1-EPS; 
	  // linear interpolation of Mu between the 2 vertexes
	  for (i=0; i<3; i++) {
	    nrm[i] = Grad[i][0][ix[1]][ix[0]]*(1-t) +
	      Grad[i][dz][ix1[1]][ix1[0]]*t;
	    // interpolation of gradient i-component between the 2 vertexes
	    x[i] = Side[i]*(-0.5*NVxl[i] + 0.5 + (float)ix[i]);
	  }
	  x[j] += Side[j]*t;
	  EdgeIpt[j][1][ix[1]][ix[0]] = StoreVertex(x, nrm);
	  // store intersection point coordinates and put its index in edgeipt
	}
      }
    }
  }
  }
  
  return 0;
}

int TriangulateGlobal::MkGrad(int iz)
//////////////////////////////////////////////////////////////////////
// Evaluates gradients in logical cube vertexes
// RETURN VALUE: zero on success 
// INPUT ARGUMENTS:
// - iz: index of slice from volumetric data
// INPUT FROM GLOBAL VARIABLES:
// **Vertex[4]: Mu Values of logical cube vertexes in 4 slices
// OUTPUT TO GLOBAL VARIABLES:
// **Grad[3][2]: gradient in cube vertexes [component][z plane][y][x]
//
// DESCRIPTION:
// This function is called by the function "Triangulate" (or "TriangulateFile")
// as part of the marching cube algorithm.
// It makes a loop on the logical cubes of a volumetric data slice. For
// each vertex of the cube it evaluates the gradient by the finit differences
// method, as described in the article
// The gradient at each logical cube vertex is stored in the global array
// **Grad[3][2]
//////////////////////////////////////////////////////////////////////
{
  int ix[3], ix0[3], ix1[3];
  int i, j, dz0, dz1;
  float Mu0, Mu1;
  
  ix[2] = iz + 1;
  for (ix[1]=0; ix[1]<NVxl[1]; ix[1]++) { // loop on cube vertexes
  for (ix[0]=0; ix[0]<NVxl[0]; ix[0]++) {
    for (i=0; i<3; i++)
      Grad[i][0][ix[1]][ix[0]] = Grad[i][1][ix[1]][ix[0]]; // update z plane
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) ix1[j] = ix0[j] = ix[j];
      if (ix0[i] > 0) ix0[i]--;
      if (ix1[i] < (NVxl[i]-1)) ix1[i]++;
      dz0 = ix[2] - ix0[2];
      dz1 = ix1[2] - ix[2];
      Mu0 = Vertex[2-dz0][ix0[1]][ix0[0]];
      Mu1 = Vertex[2+dz1][ix1[1]][ix1[0]];
      Grad[i][1][ix[1]][ix[0]] = (Mu1 - Mu0) / (Side[i]*(ix1[i]-ix0[i])) ;
    }
  }
  }

  return 0;
}

