#include "approx_tools.h"
#include <math.h>
#include <QTextStream>
#include "funcs.h"

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
// add_nodes = {x[-1], f[x[-1]], x[n], f[x[n]]}
void method_2_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F, double ax, double bx, double ay, double by,
                   double max_f, double x2, double y2, int p){
    double *Fx = new double[nx * ny];
    double *Fy = new double[nx * ny];
    double *Fxy = new double[nx * ny];
    double *temp = new double[nx];
    double *temp_ans = new double[nx];
    double *der = new double[2];
    double *add_nodes = new double[4];
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
    delete[] temp_ans;
    delete[] temp;
    temp = new double[ny];
    temp_ans = new double[ny];
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

    delete[] temp_ans;
    delete[] temp;
    temp = new double[nx];
    temp_ans = new double[nx];

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
    delete[] temp_ans;
    delete[] temp;
    delete[] der;

    // create F

    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){
            F[i + nx * (j + 16 * 0)] = fvals[i + nx * j];
            F[i + nx * (j + 16 * 1)]  = Fy[i + nx * j];
            F[i + nx * (j + 16 * 2)]  = fvals[i + nx * (j + 1)];
            F[i + nx * (j + 16 * 3)]  = Fy[i + nx * (j + 1)];

            F[i + nx * (j + 16 * 4)]  = Fx[i + nx * j];
            F[i + nx * (j + 16 * 5)]  = Fxy[i + nx * j];
            F[i + nx * (j + 16 * 6)]  = Fx[i + nx * (j + 1)];
            F[i + nx * (j + 16 * 7)]  = Fxy[i + nx * (j + 1)];

            F[i + nx * (j + 16 * 8)] = fvals[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 9)] = Fy[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 10)]  = fvals[(i + 1) + nx * (j + 1)];
            F[i + nx * (j + 16 * 11)]  = Fy[(i + 1) + nx * (j + 1)];

            F[i + nx * (j + 16 * 12)]  = Fx[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 13)]  = Fxy[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 14)]  = Fx[(i + 1) + nx * (j + 1)];
            F[i + nx * (j + 16 * 15)]  = Fxy[(i + 1) + nx * (j + 1)];
        }
    }

    delete[] Fx;
    delete[] Fy;
    delete[] Fxy;

    // create Gamma
    temp = new double[16];
    temp_ans = new double[16];
    double* temp1 = new double[16];
    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){
            matrix_product(get_a(x_vals[i + 1] - x_vals[i]), F[i + j * nx], 4, 4, 4, temp);
            temp1 = get_a(y_vals[j + 1] - y_vals[j]);
            transpose(temp1, 4);
            matrix_product(temp, temp1, 4, 4, 4, temp_ans);
            for(int k = 0; k < 16; k++){
                if(fabs(temp_ans[k]) > EPS){
                    F[i + nx * (j + 16 * k)] = temp_ans[k];
                }else{
                    F[i + nx * (j + 16 * k)]  = 0;
                }
            }
        }
    }
    delete[] temp1;
    delete[] temp;
    delete[] temp_ans;
}

void method_1_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F){
    double *Fx = new double[nx * ny];
    double *Fy = new double[nx * ny];
    double *Fxy = new double[nx * ny];
    double *temp = new double[nx];
    double *temp_ans = new double[nx];
    double *der = new double[2];
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

    // Fy create
    delete[] temp_ans;
    delete[] temp;
    temp = new double[ny];
    temp_ans = new double[ny];
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

    delete[] temp_ans;
    delete[] temp;
    temp = new double[nx];
    temp_ans = new double[nx];

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
    delete[] temp_ans;
    delete[] temp;
    delete[] der;

    // create F

    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){
            F[i + nx * (j + 16 * 0)] = fvals[i + nx * j];
            F[i + nx * (j + 16 * 1)]  = Fy[i + nx * j];
            F[i + nx * (j + 16 * 2)]  = fvals[i + nx * (j + 1)];
            F[i + nx * (j + 16 * 3)]  = Fy[i + nx * (j + 1)];

            F[i + nx * (j + 16 * 4)]  = Fx[i + nx * j];
            F[i + nx * (j + 16 * 5)]  = Fxy[i + nx * j];
            F[i + nx * (j + 16 * 6)]  = Fx[i + nx * (j + 1)];
            F[i + nx * (j + 16 * 7)]  = Fxy[i + nx * (j + 1)];

            F[i + nx * (j + 16 * 8)] = fvals[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 9)] = Fy[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 10)]  = fvals[(i + 1) + nx * (j + 1)];
            F[i + nx * (j + 16 * 11)]  = Fy[(i + 1) + nx * (j + 1)];

            F[i + nx * (j + 16 * 12)]  = Fx[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 13)]  = Fxy[(i + 1) + nx * j];
            F[i + nx * (j + 16 * 14)]  = Fx[(i + 1) + nx * (j + 1)];
            F[i + nx * (j + 16 * 15)]  = Fxy[(i + 1) + nx * (j + 1)];
        }
    }

    delete[] Fx;
    delete[] Fy;
    delete[] Fxy;

    // create Gamma
    temp = new double[16];
    temp_ans = new double[16];
    double* temp1 = new double[16];

    double *A1 = new double[16];
    double *A2 = new double[16];
    double temp_val = 0; // for matrix multiply
    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){

            get_a(A1, x_vals[i + 1] - x_vals[i]);
            get_a(A2, y_vals[j + 1] - y_vals[j]);

            for(int i = 0; col < 4; i++){
                for(int j = 0; row < 4; j++){
                    temp = 0;
                    for(int k = 0; k < 4; k++){
                        F[i + nx * (j + 16 * (k + ))]
                        temp += a[i + k * n] * b[k + m * j];
                    }
                    ans[i + n * j] = temp;
                }
            }

            matrix_product(get_a(x_vals[i + 1] - x_vals[i]), , 4, 4, 4, temp);
            temp1 = get_a(y_vals[j + 1] - y_vals[j]);
            transpose(temp1, 4);
            matrix_product(temp, temp1, 4, 4, 4, temp_ans);
            for(int k = 0; k < 16; k++){
                if(fabs(temp_ans[k]) > EPS){
                    F[i + nx * (j + 16 * k)] = temp_ans[k];
                }else{
                    F[i + nx * (j + 16 * k)] = 0;
                }
            }
        }
    }
    delete[] temp1;
    delete[] temp;
    delete[] temp_ans;
}

double method_compute(double x, double y, double *x_vals, double *y_vals, int nx, int ny, auto F){
    double ans = 0;
    for(int i = 0; i < nx - 1; i++){
        for(int j = 0; j < ny - 1; j++){
            if(x > x_vals[i] && x < x_vals[i + 1] && y > y_vals[j] && y < y_vals[j + 1]){
                double valx = x - x_vals[i];
                double valy = y - y_vals[i];
                for(int k = 0; k < 4; k++){
                    for(int l = 0; l < 4; l++){
                        ans += pow(valx, k) * pow(valy, l) * F[i + nx * (j + 16 * (k + l * 4))];
                    }
                }
            }
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

    ans[0 + 4 * 2] = -3 / (h * h);
    ans[1 + 4 * 2] = -2 / (h);
    ans[2 + 4 * 2] = 3 / (h * h);
    ans[3 + 4 * 2] = -1 / (h);

    ans[0 + 4 * 3] = 2 / (h * h * h);
    ans[1 + 4 * 3] = 1 / (h * h);
    ans[2 + 4 * 3] = -2 / (h * h * h);
    ans[3 + 4 * 3] = 1 / (h * h);
}

// a * b, a is n1 rows, m1 columns matrix, b is n2 rows, m2 columns matrix. assuming m1 == n2
void matrix_product(double *a, double *b, int n, int m, int p, double* ans){
    double temp = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < p; j++){
            temp = 0;
            for(int k = 0; k < m; k++){
                temp += a[i + k * n] * b[k + m * j];
            }
            ans[i + n * j] = temp;
        }
    }
}


// transpose NxN matrix
void transpose(double *a, int n){
    double temp = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < i; j++){
            temp = a[i + n * j];
            a[i + n * j] = a[j + n * i];
            a[i + n * j] = temp;
        }
    }
}
