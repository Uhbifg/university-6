#ifndef MPI_MATRIX_NORM_H
#define MPI_MATRIX_NORM_H
double norm(double *array, double *inverse, int n, int *vec, int shift, int rank, int total_size, double *row_buffer, double *column_buffer);
#endif