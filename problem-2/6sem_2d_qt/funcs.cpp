#include "funcs.h"
#include <QTextStream>
#include <math.h>
#include "approximate_tools.h"

#define EPS 1.e-6

double residual_1_0(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p){
    return fabs(f_0(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p) - eval_appr(x, a_border, b_border, n, a, nodes));
}

double residual_1_1(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p){
    return fabs(f_1(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p) - eval_appr(x, a_border, b_border, n, a, nodes));
}

double residual_1_2(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p){
    return fabs(f_2(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p) - eval_appr(x, a_border, b_border, n, a, nodes));
}

double residual_1_3(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p){
    return fabs(f_3(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p) - eval_appr(x, a_border, b_border, n, a, nodes));
}

double residual_1_4(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p){
    return fabs(f_4(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p) - eval_appr(x, a_border, b_border, n, a, nodes));
}

double residual_1_5(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p){
    return fabs(f_5(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p) - eval_appr(x, a_border, b_border, n, a, nodes));
}

double residual_1_6(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p){
    return fabs(f_6(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p) - eval_appr(x, a_border, b_border, n, a, nodes));
}


static double f_0(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p) {
    if(fabs(x - x2) < EPS){
        return 1 + p * 0.01 * max_f;
    }
    return 1;
}
static double f_0_der(double x) {
    return 0;
}

static double f_1(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p) {
    if(fabs(x - x2) < EPS){
        return x + p * 0.01 * max_f;
    }
    return x;
}

static double f_1_der(double x) {
    return x;
}

static double f_2(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p) {
    if(fabs(x - x2) < EPS){
        return x * x + p * 0.01 * max_f;
    }
    return x * x;
}

static double f_2_der(double x) {
    return 2 * x;
}


static double f_3(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p) {
    if(fabs(x - x2) < EPS){
        return x * x * x + p * 0.01 * max_f;
    }
    return x * x * x;
}

static double f_3_der(double x) {
    return 3 * x * x;
}


static double f_4(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p) {
    if(fabs(x - x2) < EPS){
        return x * x * x * x + p * 0.01 * max_f;
    }
    return x * x * x * x;
}

static double f_4_der(double x) {
    return 4 * x * x * x;
}


static double f_5(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p) {
    if(fabs(x - x2) < EPS){
        return exp(x) + p * 0.01 * max_f;
    }
    return exp(x);
}

static double f_5_der(double x) {
    return exp(x);
}

static double f_6(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p) {
    if(fabs(x - x2) < EPS){
        return 1 / (25 * x * x + 1) + p * 0.01 * max_f;
    }
    return 1 / (25 * x * x + 1);
}

static double f_6_der(double x) {
    return -50 * x / ((25 *
                       x * x + 1) * (25 *
                                     x * x + 1));
}

//double a_border, double b_border, int n, double *a, double *nodes
double calc_der(int func_id, double x) {
    switch (func_id) {
        case 0:
            return f_0_der(x);
        case 1:
            return f_1_der(x);
        case 2:
            return f_2_der(x);
        case 3:
            return f_3_der(x);
        case 4:
            return f_4_der(x);
        case 5:
            return f_5_der(x);
        case 6:
            return f_6_der(x);
        default:
            QTextStream out(stdout);
            out << QString("wrong func_id");
            return 0;
    }
}

double calc_f(int func_id, double x, double *nodes, double x2, double max_f, int p) {
    switch (func_id) {
        case 0:
            return f_0(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p);
        case 1:
            return f_1(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p);
        case 2:
            return f_2(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p);
        case 3:
            return f_3(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p);
        case 4:
            return f_4(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p);
        case 5:
            return f_5(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p);
        case 6:
            return f_6(x, 0, 0, 0, nullptr, nullptr, x2, max_f, p);
        default:
            QTextStream out(stdout);
            out << QString("wrong func_id");
            return 0;
    }
}

