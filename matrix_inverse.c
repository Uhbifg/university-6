#include "matrix_inverse.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include "matrix_print.h"
/* get process num by column */
int col2process(int col, int shift){
    int process = -1;
    while(col >= 0){
        col -= shift;
        process += 1;
    }
    return process;
}

int matrix_inverse(double *array, int n, double *inverse, int *vec, int shift, int rank, int total_size, double *row_buffer, double *column_buffer) {
    double eps = 0.00000001;
    if(rank == 0){
        for(int i = 0; i < n; i++){
            vec[i] = i;
        }
    }
    /* inverse init */
    for(int i = 0; i < n; i++){
        for(int j = 0; j < shift; j++){
            if(i == j + rank * shift){
                inverse[j + i * shift] = 1;
            }else{
                inverse[j + i * shift] = 0;
            }
        }
    }

#if defined DEBUG
    printf("\n inverse matrix init: \n");
    matrix_print(inverse, n, n, 1, vec, shift, rank, total_size, row_buffer);
#endif

    for(int i = 0; i < n; i++){
        int max_col = 0;
        double max_element = 0;
        int proc = 0;
        int temp = 0;
        MPI_Gather(array + i * shift, shift, MPI_DOUBLE, row_buffer, shift, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /* move i'th row to row_buffer */
        if(rank == 0){
            max_element = array[0 + shift*i];
            for(int j = 0; j < n; j++){
                if(row_buffer[j] > max_element){
                    max_col = j;
                    max_element = row_buffer[j];
                }
            }
            if(fabs(max_element) < eps){
                printf("matrix has det 0, sorry");
                return -1;
            }
            temp = vec[i];
            vec[i] = max_col;
            vec[max_col] = temp;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&max_element, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&max_col, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        proc = col2process(max_col, shift);
        if(rank == proc){
            for(int j = 0; j < n; j++){
                column_buffer[j] = array[max_col % shift + j * shift];
            }
        }
        MPI_Bcast(column_buffer, n, MPI_DOUBLE, proc, MPI_COMM_WORLD);

        for(int j = 0; j < n; j++){
            if(i != j){
                for (int k = 0; k < shift; k++) {
                    array[k + shift * j] -=  column_buffer[j] * array[k + shift * i] / max_element;
                    inverse[k + shift * j] -=  column_buffer[j] * inverse[k + shift * i] / max_element;
                }
            }
        }
        for(int k = 0; k < shift; k++){
            array[k + shift * i] /= max_element;
            inverse[k + shift * i] /= max_element;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return 0;
}
