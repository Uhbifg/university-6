#include "norm.h"
#include <math.h>
#include <mpi.h>
#include "matrix_inverse.h"
#include "stdio.h"
double norm(double *array, double *inverse, int n, int *vec, int shift, int rank, int total_size, double *row_buffer, double *column_buffer){
    double ans = 0;
    double global_ans = 0;
    for(int i = 0; i < n; i++){
        MPI_Gather(inverse + vec[n + i] * shift, shift, MPI_DOUBLE, column_buffer, shift, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(column_buffer, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if defined DEBUG
if(rank == 0){
for(int i = 0; i < n; i++){
printf("%f ", column_buffer[i]);
}
}
#endif
        for(int j = 0; j < shift; j++){
            for(int k = 0; k < n; k++){
                ans += column_buffer[k] * array[j + k * shift];
            }
        }
    }
    MPI_Reduce(&ans, &global_ans, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_ans, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    global_ans -= n;
#if defined DEBUG
printf("global ans=%f \n", global_ans);
#endif

	
    return sqrt(global_ans);
}
