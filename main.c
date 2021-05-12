#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>  // for strtol
#include <errno.h>   // for errno


#include "matrix_init.h"
#include "matrix_print.h"
//#include "matrix_inverse.h"
//#include "norm.h"


int main(int argc, char **argv) {
    int m = 0, n = 0, k = 0, total_size = 0, rank = 0;
    char *filename = NULL;
    double* mat = NULL;
    double* inverse = NULL;
    int* vec = NULL;
    //double tv1 = 0, tv2 = 0;
    double *row_buffer = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_size);

    /* read input arguments */
    if (argc == 4) {
        char *p1, *p2, *p3;
        errno = 0;

        n = strtol(argv[1], &p1, 10);
        m = strtol(argv[2], &p2, 10);
        k = strtol(argv[3], &p3, 10);

        if (errno != 0 || *p1 != '\0' || *p2 != '\0' || *p3 != '\0') {
            if(rank == 0){
                printf("Invalid argument format \n");
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else if (argc == 5) {
        char *p1, *p2, *p3;
        errno = 0;

        n = strtol(argv[1], &p1, 10);
        m = strtol(argv[2], &p2, 10);
        k = strtol(argv[3], &p3, 10);
        filename = argv[4];

        if (errno != 0 || *p1 != '\0' || *p2 != '\0' || *p3 != '\0') {
            if(rank == 0){
                printf("Invalid argument format \n");
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else {
        if(rank == 0){
            printf("Invalid argument format \n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }


    if(k < 0 || k > 4){
        if(rank == 0){
            printf("Invalid argument format \n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* end arguments reading */
    int shift = (n + total_size - 1) / total_size;
    int start_col = shift * rank;
    int end_col = start_col + shift;

    /* create matrix */
    mat = (double*)malloc(n * shift * sizeof(double));
    inverse = (double*)malloc(n * shift * sizeof(double));

    for(int i = 0; i < shift * n; i++){
        mat[i] = 0;
        inverse[i] = 0;
    }
    if (rank == total_size - 1){
        end_col = n;
    }

    if(rank == 0){
        vec = (int*)malloc(2 * n * sizeof(int));
        row_buffer = (double*)malloc(shift * total_size * sizeof(double));
        for(int i = 0; i < shift * total_size; i++){
            row_buffer[i] = 0;
        }
    }

    if(matrix_init(mat, n, k, filename, rank, total_size, start_col, end_col, row_buffer, shift) != 0){
        if(rank == 0){
            printf("Matrix init error. \n");
            free(vec);
            free(row_buffer);
        }
        free(mat);
        free(inverse);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* print matrix */
    if(rank == 0){
        printf("Matrix: \n");
    }
    matrix_print(mat, n, m, 1, vec, shift, rank, total_size, row_buffer);

    MPI_Barrier(MPI_COMM_WORLD);


  

    #if defined DEBUG
    #endif

    free(mat);
    free(inverse);
    if(rank == 0){
        free(vec);
        free(row_buffer);
    }
    MPI_Finalize();
    return 0;
}
