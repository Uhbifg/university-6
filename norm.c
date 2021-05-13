#include "norm.h"
#include <math.h>
#include <mpi.h>


double norm(double *array, double *inverse, int n, int *vec, int shift, int rank, int total_size, double *row_buffer){
    double ans = 0;
    double element = 0;
    return 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            element = 0;
            MPI_Gather(array + vec[i] * shift, shift, MPI_DOUBLE, row_buffer, shift, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            for(int k = 0; k < n; k++){
                element += array[i + k * n] * inverse[vec[n + k] + n * j];
            }
            if(i == j){
                element -= 1;
            }
            ans += element * element;
        }
    }
    return sqrt(ans);
}
