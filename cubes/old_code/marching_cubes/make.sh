gcc -Wall -O4 -o triangulate TriangulateMain.c Triangulate.c TriangulateGlobal.c FindVertexes.c OutputData.c  PolygonizeCube.c matrix.c -lm

gcc -Wall -O4 -o triangulate_test TriangulateTest.c Triangulate.c TriangulateGlobal.c FindVertexes.c OutputData.c  PolygonizeCube.c matrix.c -lm

gcc -Wall  -O4 -o trvtx2vrml trvtx2vrml_new.c -lm 