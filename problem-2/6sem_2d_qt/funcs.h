#ifndef INC_6SEM_2D_QT_FUNCS_H
#define INC_6SEM_2D_QT_FUNCS_H
static double f_0(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
static double f_0_der(double x);

static double f_1(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
static double f_1_der(double x);

static double f_2(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
static double f_2_der(double x);

static double f_3(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
static double f_3_der(double x);

static double f_4(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
static double f_4_der(double x);

static double f_5(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
static double f_5_der(double x);

static double f_6(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
static double f_6_der(double x);

double calc_der(int func_id, double x);
double calc_f(int func_id, double x, double *nodes, double x2, double max_f, int p);
#endif //INC_6SEM_2D_QT_FUNCS_H
