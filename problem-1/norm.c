#include "norm.h"
#include <math.h>
#include <mpi.h>
#include "matrix_inverse.h"

double norm(double *array, double *inverse, int n, int *vec, int shift, int rank, int total_size, double *row_buffer,
            double *column_buffer) {
    double element = 0;
    double global_ans = 0;
    double local_ans;
    for (int i = 0; i < n; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(inverse + vec[i] * shift, shift, MPI_DOUBLE, column_buffer, shift, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(column_buffer, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int j = 0; j < shift; j++) {
            element = 0;
            if (j + shift * rank >= n) {
                break;
            }
            for (int k = 0; k < n; k++) {
                element += column_buffer[k] * array[j + k * shift];
            }
            if (i == j + shift * rank) {
                element -= 1;
            }
            local_ans += element * element;
        }
    }
    MPI_Reduce(&local_ans, &global_ans, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_ans, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return sqrt(global_ans);
}
