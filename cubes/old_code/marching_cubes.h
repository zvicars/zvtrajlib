typedef struct {
   double x,y,z;
} XYZ;
typedef struct {
   XYZ p[3];
} TRIANGLE;
typedef struct {
   XYZ p[8];
   double val[8];
} GRIDCELL;

/*
   Polygonise a tetrahedron given its vertices within a cube
   This is an alternative algorithm to polygonisegrid.
   It results in a smoother surface but more triangular facets.

                      + 0
                     /|\
                    / | \
                   /  |  \
                  /   |   \
                 /    |    \
                /     |     \
               +-------------+ 1
              3 \     |     /
                 \    |    /
                  \   |   /
                   \  |  /
                    \ | /
                     \|/
                      + 2

   It's main purpose is still to polygonise a gridded dataset and
   would normally be called 6 times, one for each tetrahedron making
   up the grid cell.
   Given the grid labelling as in PolygniseGrid one would call
      PolygoniseTri(grid,iso,triangles,0,2,3,7);
      PolygoniseTri(grid,iso,triangles,0,2,6,7);
      PolygoniseTri(grid,iso,triangles,0,4,6,7);
      PolygoniseTri(grid,iso,triangles,0,6,1,2);
      PolygoniseTri(grid,iso,triangles,0,6,1,4);
      PolygoniseTri(grid,iso,triangles,5,6,1,4);
*/
int PolygoniseTri(GRIDCELL g,double iso,
   TRIANGLE *tri,int v0,int v1,int v2,int v3)
{
   int ntri = 0;
   int triindex;

   /*
      Determine which of the 16 cases we have given which vertices
      are above or below the isosurface
   */
   triindex = 0;
   if (g.val[v0] < iso) triindex |= 1;
   if (g.val[v1] < iso) triindex |= 2;
   if (g.val[v2] < iso) triindex |= 4;
   if (g.val[v3] < iso) triindex |= 8;

   /* Form the vertices of the triangles for each case */
   switch (triindex) {
   case 0x00:
   case 0x0F:
      break;
   case 0x0E:
   case 0x01:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      ntri++;
      break;
   case 0x0D:
   case 0x02:
      tri[0].p[0] = VertexInterp(iso,g.p[v1],g.p[v0],g.val[v1],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      ntri++;
      break;
   case 0x0C:
   case 0x03:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      tri[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      ntri++;
      tri[1].p[0] = tri[0].p[2];
      tri[1].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      tri[1].p[2] = tri[0].p[1];
      ntri++;
      break;
   case 0x0B:
   case 0x04:
      tri[0].p[0] = VertexInterp(iso,g.p[v2],g.p[v0],g.val[v2],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v1],g.val[v2],g.val[v1]);
      tri[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      ntri++;
      break;
   case 0x0A:
   case 0x05:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      ntri++;
      tri[1].p[0] = tri[0].p[0];
      tri[1].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      tri[1].p[2] = tri[0].p[1];
      ntri++;
      break;
   case 0x09:
   case 0x06:
      tri[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      tri[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      tri[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      ntri++;
      tri[1].p[0] = tri[0].p[0];
      tri[1].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      tri[1].p[2] = tri[0].p[2];
      ntri++;
      break;
   case 0x07:
   case 0x08:
      tri[0].p[0] = VertexInterp(iso,g.p[v3],g.p[v0],g.val[v3],g.val[v0]);
      tri[0].p[1] = VertexInterp(iso,g.p[v3],g.p[v2],g.val[v3],g.val[v2]);
      tri[0].p[2] = VertexInterp(iso,g.p[v3],g.p[v1],g.val[v3],g.val[v1]);
      ntri++;
      break;
   }

   return(ntri);
}
