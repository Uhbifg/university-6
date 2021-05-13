#include "matrix_inverse.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include "matrix_print.h"
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
    //int max_col = 0;
    double max_element = 0;
    for(int i = 0; i < n; i++){
MPI_Gather(array + i * shift, shift, MPI_DOUBLE, row_buffer, shift, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /* move i'th row to row_buffer */
        if(rank == 0){
            for(int j = 0; j < n; j++){
                if(row_buffer[j] > max_element){
                    //max_col = j;
                    max_element = row_buffer[j];

                }
printf("%f ", row_buffer[j]);
            }
            if(fabs(max_element) < eps){
                printf("matrix has det 0, sorry");
                return -1;
            }
	
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&max_element, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        printf("max_el %f", max_element);
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
        MPI_Barrier(MPI_COMM_WORLD);
#if defined DEBUG
        if(rank == 0){
            printf("step: %i \n matrix: \n", i);
        }

        matrix_print(array, n, n, 1, vec, shift, rank, total_size, row_buffer);
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    return 0;
}
