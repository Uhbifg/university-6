#include "matrix_print.h"
#include <stdio.h>
#include <mpi.h>


void matrix_print(double *array, int n, int m, int flag, int *vec, int shift, int rank, int total_size, double *row_buffer){
    for(int i = 0; i < n; i++){
        /* end if m rows reached */
        if (i == m){
            break;
        }
        if(flag == 1){
#if defined DEBUG_print

            if(rank == 0){
		printf("\n row_buff before: \n");
		    for(int i = 0; i < total_size * shift; i++){
		        printf("%10.3e ", row_buffer[i]);
		    }
		    printf("\n");
		}
#endif
            MPI_Gather(array + i * shift, shift, MPI_DOUBLE, row_buffer, shift, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined DEBUG_print

            if(rank == 0){
			printf("\n row_buff after: \n");
		    for(int i = 0; i < total_size * shift; i++){
		        printf("%10.3e ", row_buffer[i]);
		    }
		    printf("\n %i \n", shift);
		}
#endif
        }else{
            MPI_Gather(array + vec[i] * shift, shift, MPI_DOUBLE, row_buffer, shift, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if(rank == 0){
            for(int j = 0; j < n; j++){
                if(j == m){
                    break;
                }
                if(j != n - 1){
                    printf("%10.3e ", row_buffer[j]);
                }else {
                    printf("%10.3e", row_buffer[j]);
                }
            }
            printf("\n");
        }

    }
    if(rank == 0){
        printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
