#include "matrix_init.h"

#include "matrix_init.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

double f(int k, int n, int i, int j) {
    switch (k) {
        case 1:
            if (i > j) {
                return n - i;
            } else {
                return n - j;
            }
        case 2:
            if (i > j) {
                return i + 1;
            } else {
                return j + 1;
            }
        case 3:
            if (i > j) {
                return i - j;
            } else {
                return j - i;
            }
        case 4:
            return 1 / (i + j + 1.0);
    }
    return 0;
}


int matrix_init(double *array, int n, int k, char *filename, int rank, int total_size, int start_col, int end_col, double *row_buffer, int shift) {
    FILE *inp;
    if (k == 0) {
        if(rank == 0){
            inp = fopen(filename, "r");
            if (inp == NULL) {
                printf("Cannot open file.\n");
                return -1;
            }
        }
        for (int i = 0; i < n; i++) {
            if(rank == 0) {
                for (int j = 0; j < n; j++) {
                    if (fscanf(inp, "%lf,", &row_buffer[j]) != 1) {
                        printf("File read error. \n");
                        fclose(inp);
                        return -1;
                    }
                    if (j != n - 1) {
                        if (fgetc(inp) != ' ') {
                            printf("Invalid matrix format. \n");
                            fclose(inp);
                            return -1;
                        }
                    }
                }
                if (fgetc(inp) != '\n') {
                    printf("Invalid matrix format. \n");
                    fclose(inp);
                    return -1;
                }
            }
            MPI_Scatter(row_buffer, shift, MPI_DOUBLE, array + i * shift, shift, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
        }
        fclose(inp);
    } else {
        for (int i = 0; i < n; i++) {
            for (int j = start_col; j < end_col; j++) {
                array[j - start_col + shift * i] = f(k, n, i, j);
            }
        }
    }
    return 0;
}
