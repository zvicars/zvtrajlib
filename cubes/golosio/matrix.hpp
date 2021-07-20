#ifndef MATRIXH
#define MATRIXH

// allocate a float matrix with dimensions nx, ny
float **float_matrix(long nx, long ny);

// allocate a int matrix with dimensions nx, ny
int **int_matrix(long nx, long ny);

// free a float matrix allocated by float_matrix()
int free_float_matrix(float **m);

// free a int matrix allocated by int_matrix()
int free_int_matrix(int **m);

#endif
