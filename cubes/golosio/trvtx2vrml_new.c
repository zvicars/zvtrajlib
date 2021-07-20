/****************************************************/
/*                                                  */
/*       Virtual Reality Modelling Language         */
/*               CONVERTER                          */
/*                                                  */
/****************************************************/
/*        To oonvert raw files in vrml files        */
/****************************************************/


#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
   FILE *tr_fp, *vtx_fp, *vrml_fp;
   int n_tr, n_vtx;
   float vtx[6];
   int iv[3];
   int i;
   
   if (argc != 4) {
      printf("Usage: %s tr-input-file vtx-input-file vrml-output-file\n", argv[0]);
      printf("Use extention .wrl for  vrml-output-file\n");
      exit(0);
   }
   if ((vrml_fp=fopen(argv[3], "w"))==NULL) {
      printf("Cannot open output file\n");
      exit(EXIT_FAILURE);
   }
   if ((vtx_fp=fopen(argv[2], "rb"))==NULL) {
      printf("Cannot open vertex file\n");
      exit(EXIT_FAILURE);
   }
   fread(&n_vtx, sizeof(int), 1, vtx_fp);
   //fprintf(vrml_fp, "%d\n", n_vtx);

   /*template for vrml files*/
   fprintf(vrml_fp,"#VRML V2.0 utf8\n");  
   fprintf(vrml_fp,"Shape {\n"); 
   fprintf(vrml_fp,"appearance Appearance {\n"); 
   fprintf(vrml_fp,"material Material { emissiveColor .8 0.2 0 }\n}\n");
   fprintf(vrml_fp,"geometry IndexedFaceSet{\n");
   fprintf(vrml_fp,"coord Coordinate {\n");
   fprintf(vrml_fp,"point[\n");

 for (i=0; i<n_vtx; i++) {
      fread(vtx, sizeof(float), 6, vtx_fp);
      fprintf(vrml_fp, " %.3e %.3e %.3e", vtx[0], vtx[1], vtx[2]);
      if (i<n_vtx-1) fprintf(vrml_fp, ",\n");
   }
   
   fclose(vtx_fp);

   fprintf(vrml_fp,"]\n}\n");
   fprintf(vrml_fp,"coordIndex [\n");

   if ((tr_fp=fopen(argv[1], "rb"))==NULL) {
      printf("Cannot open triangle file\n");
      exit(EXIT_FAILURE);
   }
   fread(&n_tr, sizeof(int), 1, tr_fp);
   
   for (i=0; i<n_tr; i++) {
      fread(iv, sizeof(int), 3, tr_fp);
      fprintf(vrml_fp, " %d, %d, %d, -1\n", iv[0], iv[1], iv[2]);
   }
   fclose(tr_fp);
   
   /*ADD gouraud shading correction*/
   fprintf(vrml_fp,"]\ncreaseAngle 0.1\n");
   fprintf(vrml_fp,"}\n}\n");

   fclose(vrml_fp);
   
   return 0;
}
