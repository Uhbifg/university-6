#include "matrix_inverse.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
/* get process num by column */
int col2process(int col, int shift, int n){
    int process = -1;
    while(col > 0){
        col -= shift;
        process += 1;
    }
    return process;
}

int matrix_inverse(double *array, int n, double *inverse, int *vec, int shift, int rank, int total_size, double *row_buffer) {
    double eps = 0.00000001;
    int max_col = 0;
    double max_element = 0;
    for(int i = 0; i < n; i++){
        /* move i'th row to row_buffer */
        if(rank == 0){
            MPI_Scatter(row_buffer, shift, MPI_DOUBLE, array + i * shift, shift, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
            for(int j = 0; j < n; j++){
                if(row_buffer[j] > row_buffer[max_col]){
                    max_col = j;
                    max_element = row_buffer[max_col];
                    if(fabs(max_element) < eps){
                        printf("matrix has det 0, sorry");
                        return -1;
                    }
                }
            }
            MPI_Bcast(&max_element, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        for(int j = 0; j < n; j++){
            for(int k = 0; k < shift; k++){
                if(j != i){
                    array[k + shift * j] -= array[k + shift * i] / max_element;
                }
            }
        }

        for(int k = 0; k < shift; k++){
            array[k + shift * i] /= max_element;
        }
    }
    return 0;
}
