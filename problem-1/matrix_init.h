//
// Created by admin on 12.05.2021.
//

#ifndef MPI_MATRIX_MATRIX_INIT_H
#define MPI_MATRIX_MATRIX_INIT_H
int matrix_init(double *array, int n, int k, char *filename, int rank, int total_size, int start_col, int end_col, double *row_buffer, int shift);
double f(int k, int n, int i, int j);

#endif //MPI_MATRIX_MATRIX_INIT_H
