#include "approx_tools.h"
#include <math.h>
#include <QTextStream>
#include "funcs.h"
#define MAX(a, b) ((a) < (b) ? (b) : (a))

#define EPS 1e-6
int sign(double x) {
    if (x > 0) {
        return 1;
    } else if (x < 0) {
        return -1;
    }
    return 0;
}

double min(double x1, double x2) {
    if (x1 < x2) {
        return x1;
    }
    return x2;
}

double max(double x1, double x2) {
    if (x1 < x2) {
        return x2;
    }
    return x1;
}

int method_init_1(int n, double *x, double *f_vals, double *additional_space, double *der) {
    double w1 = 0, w2 = 0;
    double eps = 1e-6;
    for (int i = 2; i <= n - 3; i++) {
        w1 = w_j(f_vals, x, i + 1);
        w2 = w_j(f_vals, x, i - 1);
        if (w1 < eps && w2 < eps) {
            additional_space[i] =
                    ((x[i + 1] - x[i]) * dd(f_vals, x, i - 1, i) + (x[i] - x[i - 1]) * dd(f_vals, x, i, i + 1)) /
                    (x[i + 1] - x[i - 1]);
        } else {
            additional_space[i] = (w1 * dd(f_vals, x, i - 1, i) + w2 * dd(f_vals, x, i, i + 1)) / (w1 + w2);
        }
    }

    additional_space[0] = der[0]; // f'(x_0)
    additional_space[n - 1] = der[1]; // f'(x_{n-1})

    additional_space[1] = 3 * dd(f_vals, x, 0, 1) - 2 * additional_space[0];
    additional_space[n - 2] = 3 * dd(f_vals, x, n - 1, n - 2) - 2 * additional_space[n - 1];

    return 0;
}


double w_j(double *f_vals, double *x, int j) {
    return fabs(dd(f_vals, x, j, j + 1) - dd(f_vals, x, j - 1, j));
}

// divided differences f(x_i; x_j) = (f(x_j) - f(x_i)) / (x_j - x_i)
double dd(double *f_vals, double *x, int i, int j) {
    return (f_vals[j] - f_vals[i]) / (x[j] - x[i]);
}




// add_nodes = {x[-1], f[x[-1]], x[n], f[x[n]]}
int method_init_2(int n, double *x, double *f_vals, double *additional_space, double *add_nodes) {
    double dd1 = 0, dd2 = 0;
    for (int i = 1; i <= n - 2; i++) {
        if (i == 1) {
            dd1 = dd(f_vals, x, i - 1, i);
            dd2 = dd(f_vals, x, i, i + 1);
        } else {
            dd1 = dd2;
            dd2 = dd(f_vals, x, i, i + 1);
        }

        if (sign(dd1) * sign(dd2) == 1) {
            additional_space[i] = sign(dd2) * min(fabs(dd1), fabs(dd2));
        } else {
            additional_space[i] = 0;
        }
    }


    dd1 = (f_vals[0] - add_nodes[1]) / (x[0] - add_nodes[0]);
    dd2 = dd(f_vals, x, 0, 1);
    if (sign(dd1) * sign(dd2) == 1) {
        additional_space[0] = sign(dd2) * min(fabs(dd1), fabs(dd2));
    } else {
        additional_space[0] = 0;
    }

    dd1 = dd(f_vals, x, n - 2, n - 1);
    dd2 = (add_nodes[3] - f_vals[n - 1]) / (add_nodes[2] - x[n - 1]);
    if (sign(dd1) * sign(dd2) == 1) {
        additional_space[n - 1] = sign(dd2) * min(fabs(dd1), fabs(dd2));
    } else {
        additional_space[n - 1] = 0;
    }

    return 0;
}
void method_1_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F){
    double *Fx = (double*) malloc(sizeof(double) * nx * ny);
    double *Fy = (double*) malloc(sizeof(double) * nx * ny);
    double *Fxy = (double*) malloc(sizeof(double) * nx * ny);
    double *temp = (double*) malloc(sizeof(double) * MAX(nx, ny));
    double *temp_ans = (double*) malloc(sizeof(double) * MAX(nx, ny));
    double *der = (double*) malloc(sizeof(double) * 2);;
    // Fx create
    for(int i = 0; i < ny; i++){
        for(int j = 0; j < nx; j++){
            temp[j] = fvals[j + nx * i];
        }
        der[0] = calc_der_x(x_vals[0], y_vals[i], func_id);
        der[1] = calc_der_x(x_vals[nx - 1], y_vals[i], func_id);
        method_init_1(nx, x_vals, temp, temp_ans, der);
        for(int j = 0; j < nx; j++){
            Fx[j + nx * i] = temp_ans[j];
        }
    }

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            temp[j] = fvals[i + nx * j];
        }
        der[0] = calc_der_y(x_vals[i], y_vals[0], func_id);
        der[1] = calc_der_x(x_vals[i], y_vals[ny - 1], func_id);
        method_init_1(ny, y_vals, temp, temp_ans, der);
        for(int j = 0; j < ny; j++){
            Fy[i + nx * j] = temp_ans[j];
        }
    }

    // create Fxy
    for(int i = 0; i < ny; i++){
        for(int j = 0; j < nx; j++){
            temp[j] = Fy[j + nx * i];
        }
        // i don't know what else
        der[0] = 0;
        der[1] = 0;
        method_init_1(nx, x_vals, temp, temp_ans, der);
        for(int j = 0; j < nx; j++){
            Fxy[j + nx * i] = temp_ans[j];
        }
    }

    // create F -> gamma

    create_gamma(nx, ny, F, Fx, Fy, Fxy, fvals, y_vals, x_vals, func_id);

    free(temp_ans);
    free(temp);
    free(der);
    free(Fx);
    free(Fy);
    free(Fxy);
}


// add_nodes = {x[-1], f[x[-1]], x[n], f[x[n]]}
void method_2_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F, double ax, double bx, double ay, double by,
                   double max_f, double x2, double y2, int p){
    double *Fx = (double*) malloc(sizeof(double) * nx * ny);
    double *Fy = (double*) malloc(sizeof(double) * nx * ny);
    double *Fxy = (double*) malloc(sizeof(double) * nx * ny);
    double *temp = (double*) malloc(sizeof(double) * MAX(nx, ny));
    double *temp_ans = (double*) malloc(sizeof(double) * MAX(nx, ny));
    double *add_nodes = (double*) malloc(sizeof(double) * 4);
    // Fx create
    for(int i = 0; i < ny; i++){
        for(int j = 0; j < nx; j++){
            temp[j] = fvals[j + nx * i];
        }
        add_nodes[0] = (bx - ax) / (nx - 1) * (-1) + ax;
        add_nodes[2] = (bx - ax) / (nx - 1) * nx + ax;
        add_nodes[1] = calc_func(add_nodes[0], y_vals[i], func_id, max_f, p, x2, y2);
        add_nodes[3] = calc_func(add_nodes[2], y_vals[i], func_id, max_f, p, x2, y2);
        method_init_2(nx, x_vals, temp, temp_ans, add_nodes);
        for(int j = 0; j < nx; j++){
            Fx[j + nx * i] = temp_ans[j];
        }
    }

    // Fy create
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            temp[j] = fvals[i + nx * j];
        }
        add_nodes[0] = (by - ay) / (ny - 1) * (-1) + ay;
        add_nodes[2] = (by - ay) / (ny - 1) * ny + ay;
        add_nodes[1] = calc_func(x_vals[i], add_nodes[0] , func_id, max_f, p, x2, y2);
        add_nodes[3] = calc_func(x_vals[i], add_nodes[2], func_id, max_f, p, x2, y2);
        method_init_2(ny, y_vals, temp, temp_ans, add_nodes);

        for(int j = 0; j < ny; j++){
            Fy[i + nx * j] = temp_ans[j];
        }
    }


    for(int i = 0; i < ny; i++){
        for(int j = 0; j < nx; j++){
            temp[j] = Fy[j + nx * i];
        }
        // i don't know what else
        add_nodes[0] = 0;
        add_nodes[2] = 0;
        add_nodes[1] = 0;
        add_nodes[3] = 0;
        method_init_2(nx, x_vals, temp, temp_ans, add_nodes);
        for(int j = 0; j < nx; j++){
            Fxy[j + nx * i] = temp_ans[j];
        }
    }
    // create F -> gamma
    create_gamma(nx, ny, F, Fx, Fy, Fxy, fvals, y_vals, x_vals, func_id);
    /*
    free(temp_ans);
    free(temp);
    free(add_nodes);
    free(Fx);
    free(Fy);
    free(Fxy);
*/
     }

int create_gamma(int nx, int ny, double*F, double *Fx, double *Fy, double *Fxy, double *fvals, double *y_vals, double *x_vals, int func_id){
    double temp_val = 0; // for matrix multiply

    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){
            F[i * 16 * ny + j * 16 + 0] = fvals[i + nx * j];
            F[i * 16 * ny + j * 16 + 1]  = Fy[i + nx * j];
            F[i * 16 * ny + j * 16 + 2]  = fvals[i + nx * (j + 1)];
            F[i * 16 * ny + j * 16 + 3]  = Fy[i + nx * (j + 1)];

            F[i * 16 * ny + j * 16 + 4]  = Fx[i + nx * j];
            F[i * 16 * ny + j * 16 + 5]  = Fxy[i + nx * j];
            F[i * 16 * ny + j * 16 + 6]  = Fx[i + nx * (j + 1)];
            F[i * 16 * ny + j * 16 + 7]  = Fxy[i + nx * (j + 1)];

            F[i * 16 * ny + j * 16 + 8] = fvals[(i + 1) + nx * j];
            F[i * 16 * ny + j * 16 + 9] = Fy[(i + 1) + nx * j];
            F[i * 16 * ny + j * 16 + 10]  = fvals[(i + 1) + nx * (j + 1)];
            F[i * 16 * ny + j * 16 + 11]  = Fy[(i + 1) + nx * (j + 1)];

            F[i * 16 * ny + j * 16 + 12]  = Fx[(i + 1) + nx * j];
            F[i * 16 * ny + j * 16 + 13]  = Fxy[(i + 1) + nx * j];
            F[i * 16 * ny + j * 16 + 14]  = Fx[(i + 1) + nx * (j + 1)];
            F[i * 16 * ny + j * 16 + 15]  = Fxy[(i + 1) + nx * (j + 1)];
        }
    }

    double* temp = (double*) malloc(sizeof(double) * 16);
    double* temp_ans = (double*) malloc(sizeof(double) * 16);

    double *A1 = (double*) malloc(sizeof(double) * 16);
    double *A2 = (double*) malloc(sizeof(double) * 16);
    /*
    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){
            F[i * 16 * ny + j * 16 + 0] = calc_func(x_vals[i], y_vals[j], func_id, 0, 0, 0, 0);
            F[i * 16 * ny + j * 16 + 1]  = calc_der_y(x_vals[i], y_vals[j], func_id);
            F[i * 16 * ny + j * 16 + 2]  = calc_func(x_vals[i], y_vals[j + 1], func_id, 0, 0, 0, 0);
            F[i * 16 * ny + j * 16 + 3]  = calc_der_y(x_vals[i], y_vals[j + 1], func_id);

            F[i * 16 * ny + j * 16 + 4]  = calc_der_x(x_vals[i], y_vals[j], func_id);
            F[i * 16 * ny + j * 16 + 5]  = calc_der_xy(x_vals[i], y_vals[j], func_id);
            F[i * 16 * ny + j * 16 + 6]  = calc_der_x(x_vals[i], y_vals[j + 1], func_id);
            F[i * 16 * ny + j * 16 + 7]  = calc_der_xy(x_vals[i], y_vals[j + 1], func_id);

            F[i * 16 * ny + j * 16 + 8] = calc_func(x_vals[i], y_vals[j], func_id, 0, 0, 0, 0);
            F[i * 16 * ny + j * 16 + 9] = calc_der_y(x_vals[i + 1], y_vals[j], func_id);
            F[i * 16 * ny + j * 16 + 10]  = calc_func(x_vals[i + 1], y_vals[j + 1], func_id, 0, 0, 0, 0);
            F[i * 16 * ny + j * 16 + 11]  = calc_der_y(x_vals[i + 1], y_vals[j + 1], func_id);

            F[i * 16 * ny + j * 16 + 12]  = calc_der_x(x_vals[i + 1], y_vals[j], func_id);
            F[i * 16 * ny + j * 16 + 13]  = calc_der_xy(x_vals[i + 1], y_vals[j], func_id);
            F[i * 16 * ny + j * 16 + 14]  = calc_der_x(x_vals[i + 1], y_vals[j + 1], func_id);
            F[i * 16 * ny + j * 16 + 15]  = calc_der_xy(x_vals[i + 1], y_vals[j + 1], func_id);
        }
    }
     */
    // create gamma
    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){

            get_a(A1, x_vals[i + 1] - x_vals[i]);
            get_a(A2, y_vals[j + 1] - y_vals[j]);
            //transpose(A2, 4);

            for(int col = 0; col < 4; col++){ // i
                for(int row = 0; row < 4; row++){ // j
                    temp_val = 0;
                    for(int k = 0; k < 4; k++){
                        temp_val += A1[col + k * 4] * F[i * 16 * ny + 16 * j + k + 4 * row];
                    }
                    temp[col + 4 * row] = temp_val;
                }
            }

            for(int col = 0; col < 4; col++){
                for(int row = 0; row < 4; row++){
                    temp_val = 0;
                    for(int k = 0; k < 4; k++){
                        temp_val += temp[col + k * 4] * A2[row + 4 * k];
                    }
                    temp_ans[col + 4 * row] = temp_val;
                }
            }
            for(int k = 0; k < 16; k++){
                if(fabs(temp_ans[k]) < EPS){
                    F[i * 16 * ny + 16 * j + k] = 0;
                }else{
                    F[i * 16 * ny + 16 * j + k] = temp_ans[k];
                }
            }
        }
    }
    return 0;
}

double method_compute(double x, double y, double *x_vals, double *y_vals, int nx, int ny, double*  F){
    double ans = 0;
    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){
            if(x > x_vals[i] && x < x_vals[i + 1] && y > y_vals[j] && y < y_vals[j + 1]){
                double valx = x - x_vals[i];
                double valy = y - y_vals[i];
                for(int k = 0; k < 4; k++){
                    for(int l = 0; l < 4; l++){
                        ans += pow_d(valx, k) * pow_d(valy, l) * F[i * 16 * ny + 16 * j + (l + k * 4)];
                    }
                }
                return ans;
            }else if(fabs(x - x_vals[i]) < EPS && fabs(y - y_vals[j]) < EPS){
                return F[i * 16 * nx + 16 * j];
            }
        }
    }
    return ans;
}

double method_compute(double x, double y, double *x_vals, double *y_vals, int nx, int ny, double*  F, int i, int j){
    double ans = 0;
    double valx = x - x_vals[i];
    double valy = y - y_vals[i];
    for(int k = 0; k < 4; k++){
        for(int l = 0; l < 4; l++){
            ans += pow_d(valx, k) * pow_d(valy, l) * F[i * 16 * ny + 16 * j + (k + l * 4)];
        }
    }
    return ans;
}

void get_a(double *ans, double h){
    for(int i = 0; i < 16; i++){
        ans[i] = 0;
    }
    ans[0 + 4 * 0] = 1;

    ans[1 + 4 * 1] = 1;

    ans[2 + 4 * 0] = -3 / (h * h);
    ans[2 + 4 * 1] = -2 / (h);
    ans[2 + 4 * 2] = 3 / (h * h);
    ans[2 + 4 * 3] = -1 / (h);

    ans[3 + 4 * 0] = 2 / (h * h * h);
    ans[3 + 4 * 1] = 1 / (h * h);
    ans[3 + 4 * 2] = -2 / (h * h * h);
    ans[3 + 4 * 3] = 1 / (h * h);
}




double pow_d(double x, int p){
    double ans = 1;
    while(p != 0){
        p--;
        ans *= x;
    }
    return ans;
}
