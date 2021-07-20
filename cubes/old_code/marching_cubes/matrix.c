#include <stdio.h>
#include <stdlib.h>

// allocate a float matrix with dimensions nx, ny  
float **float_matrix(long nx, long ny)
{
  long i;
  float **m;
  
  // allocate pointers to rows
  m=(float **) malloc(sizeof(float*)*nx);
  if (!m) {
    printf("allocation error in float_matrix()");
    exit(EXIT_FAILURE);
  }
  // allocate rows and set pointers to them
  m[0]=(float *) malloc(sizeof(float)*nx*ny);
  if (!m[0]) {
    printf("allocation error in float_matrix()");
    exit(EXIT_FAILURE);
  }
  for(i=1; i<nx; i++) m[i]=m[i-1]+ny;
  // return pointer
  return m;
}

// allocate a int matrix with dimensions nx, ny  
int **int_matrix(long nx, long ny)
{
  long i;
  int **m;
  
  // allocate pointers to rows
  m=(int **) malloc(sizeof(int*)*nx);
  if (!m) {
    printf("allocation error in int_matrix()");
    exit(EXIT_FAILURE);
  }
  // allocate rows and set pointers to them
  m[0]=(int *) malloc(sizeof(int)*nx*ny);
  if (!m[0]) {
    printf("allocation error in int_matrix()");
    exit(EXIT_FAILURE);
  }
  for(i=1; i<nx; i++) m[i]=m[i-1]+ny;
  // return pointer
  return m;
}

// free a float matrix allocated by float_matrix()
int free_float_matrix(float **m)
{
        free(m[0]);
        free(m);

	return 0;
}

// free a int matrix allocated by int_matrix()
int free_int_matrix(int **m)
{
        free(m[0]);
        free(m);

	return 0;
}

