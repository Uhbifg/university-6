#include "approximate_tools.h"
#include <math.h>
#include <QTextStream>

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


int method_init_1(int n, double *x, double *f_vals, double *a, double *additional_space, double *der) {
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

    recover_coef_from_d(f_vals, x, n, a, additional_space);
    return 0;
}

void recover_coef_from_d(double *f_vals, double *x, int n, double *a, double *d) {
    for (int i = 0; i < n - 1; i++) {
        a[0 + i * 4] = f_vals[i];
        a[1 + i * 4] = d[i];
        a[2 + i * 4] = (3*dd(f_vals, x, i, i + 1) - 2 * d[i] - d[i + 1])/(x[i + 1] - x[i]);
        a[3 + i * 4] = (d[i] + d[i + 1] - 2 * dd(f_vals, x, i, i + 1)) / ((x[i + 1] - x[i]) * (x[i + 1] - x[i]));
    }
}

double w_j(double *f_vals, double *x, int j) {
    return fabs(dd(f_vals, x, j, j + 1) - dd(f_vals, x, j - 1, j));
}

// divided differences f(x_i; x_j) = (f(x_j) - f(x_i)) / (x_j - x_i)
double dd(double *f_vals, double *x, int i, int j) {
    return (f_vals[j] - f_vals[i]) / (x[j] - x[i]);
}




// add_nodes = {x[-1], f[x[-1]], x[n], f[x[n]]}
int method_init_2(int n, double *x, double *f_vals, double *a, double *additional_space, double *add_nodes) {
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

    recover_coef_from_d(f_vals, x, n, a, additional_space);
    return 0;
}

double eval_appr(double x, double a_border, double b_border, int n, double *a, double *nodes) {
    QTextStream out(stdout);
    for (int i = 0; i < n - 1; i++) {
        if (x >= nodes[i] && x <= nodes[i + 1]) {
            double y = x - nodes[i];
            return a[0 + i * 4] + a[1 + i * 4] * y + a[2 + i * 4] * y * y + a[3 + i * 4] * y * y * y;
        }
    }
    double y = x - nodes[n-1];
    return a[0 + (n-1) * 4] + a[1 + (n-1) * 4] * y + a[2 + (n-1) * 4] * y * y + a[3 + (n-1) * 4] * y * y * y;
}
