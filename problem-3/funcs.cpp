#ifndef PROBLEM_3_FUNCS_H
#define PROBLEM_3_FUNCS_H

#include "funcs.h"
#include "math.h"
double eps = 1e-9;

double f0(double x, double y) { return 1; }

double f1(double x, double y) { return x; }

double f2(double x, double y) { return y; }

double f3(double x, double y) { return x + y; }

double f4(double x, double y) { return sqrt(x * x + y * y); }

double f5(double x, double y) { return x * x + y * y; }

double f6(double x, double y) { return exp(x * x - y * y); }

double f7(double x, double y) { return 1 / (25 * (x * x + y * y) + 1); }

double Dxf0(double x, double y) { return 0; }

double Dxf1(double x, double y) { return 1; }

double Dxf2(double x, double y) { return 0; }

double Dxf3(double x, double y) { return 1; }

double Dxf4(double x, double y) {
    if (x * x + y * y > eps)
        return x / sqrt(x * x + y * y);
    else return 0;
}

double Dxf5(double x, double y) { return 2 * x; }

double Dxf6(double x, double y) { return 2 * x * exp(x * x - y * y); }

double Dxf7(double x, double y) { return -50 * x / (25 * (x * x + y * y) + 1) / (25 * (x * x + y * y) + 1); }

double Dyf0(double x, double y) { return 0; }

double Dyf1(double x, double y) { return 0; }

double Dyf2(double x, double y) { return 1; }

double Dyf3(double x, double y) { return 1; }

double Dyf4(double x, double y) {
    if (x * x + y * y > eps)
        return y / sqrt(x * x + y * y);
    else return 0;
}

double Dyf5(double x, double y) { return 2 * y; }

double Dyf6(double x, double y) { return -2 * y * exp(x * x - y * y); }

double Dyf7(double x, double y) { return -50 * y / (25 * (x * x + y * y) + 1) / (25 * (x * x + y * y) + 1); }

double DxDyf0(double x, double y) { return 0; }

double DxDyf1(double x, double y) { return 0; }

double DxDyf2(double x, double y) { return 0; }

double DxDyf3(double x, double y) { return 0; }

double DxDyf4(double x, double y) {
    if (x * x + y * y > eps)
        return -x * y / sqrt(x * x + y * y) / (x * x + y * y);
    else return 0;
}

double DxDyf5(double x, double y) { return 0; }

double DxDyf6(double x, double y) { return -4 * x * y * exp(x * x - y * y); }

double DxDyf7(double x, double y) {
    return 5000 * y * x / (25 * (x * x + y * y) + 1) / (25 * (x * x + y * y) + 1) / (25 * (x * x + y * y) + 1);
}

double calc_der_x(double x, double y, int func_id){
    double f = 0;
    switch (func_id) {
        case 0:
            f = Dxf0(x, y);
            break;
        case 1:
            f = Dxf1(x, y);
            break;
        case 2:
            f = Dxf2(x, y);
            break;
        case 3:
            f = Dxf3(x, y);
            break;
        case 4:
            f = Dxf4(x, y);
            break;
        case 5:
            f = Dxf5(x, y);
            break;
        case 6:
            f = Dxf6(x, y);
            break;
        case 7:
            f = Dxf7(x, y);
            break;
    }
    return f;
}


double calc_der_y(double x, double y, int func_id){
    double f = 0;
    switch (func_id) {
        case 0:
            f = Dyf0(x, y);
            break;
        case 1:
            f = Dyf1(x, y);
            break;
        case 2:
            f = Dyf2(x, y);
            break;
        case 3:
            f = Dyf3(x, y);
            break;
        case 4:
            f = Dyf4(x, y);
            break;
        case 5:
            f = Dyf5(x, y);
            break;
        case 6:
            f = Dyf6(x, y);
            break;
        case 7:
            f = Dyf7(x, y);
            break;
    }
    return f;
}

double calc_func(double x, double y, int func_id, double max_f, int p, double x2, double y2) {
    double f;
    switch (func_id) {
        case 0:
            f = f0(x, y);
            break;
        case 1:
            f = f1(x, y);
            break;
        case 2:
            f = f2(x, y);
            break;
        case 3:
            f = f3(x, y);
            break;
        case 4:
            f = f4(x, y);
            break;
        case 5:
            f = f5(x, y);
            break;
        case 6:
            f = f6(x, y);
            break;
        case 7:
            f = f7(x, y);
            break;
    }
    if (p != 0 && fabs(x - x2) < 1e-6 && fabs(y - y2) < 1e-6){
        f += p * 0.1 * max_f;
    }
    return f;
}
#endif //PROBLEM_3_FUNCS_H
