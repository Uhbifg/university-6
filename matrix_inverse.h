#ifndef MPI_MATRIX_MATRIX_INVERSE_H
#define MPI_MATRIX_MATRIX_INVERSE_H
int col2process(int col, int shift);
int matrix_inverse(double *array, int n, double *inverse, int *vec, int shift, int rank, int total_size, double *row_buffer, double *column_buffer);
#endif //MPI_MATRIX_MATRIX_INVERSE_H
