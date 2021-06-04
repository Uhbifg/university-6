#ifndef INC_6SEM_2D_QT_APPROXIMATE_TOOLS_H
#define INC_6SEM_2D_QT_APPROXIMATE_TOOLS_H

int sign(double x);
int method_init_1(int n, double *x, double *f_vals, double *a, double *additional_space, double *der);
double dd(double *f_vals, double *x, int i, int j);
void recover_coef_from_d(double *f_vals, double *x, int n, double *a, double *d);
int method_init_2(int n, double *x, double *f_vals, double *a, double *additional_space, double *add_nodes);
double w_j(double *f_vals, double *x, int j);
double eval_appr(double x, double a_border, double b_border, int n, double *a, double *nodes);
double max(double x1, double x2);
double residual_1_0(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
double residual_1_1(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
double residual_1_2(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
double residual_1_3(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
double residual_1_4(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
double residual_1_5(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);
double residual_1_6(double x, double a_border, double b_border, int n, double *a, double *nodes, double x2, double max_f, int p);

#endif //INC_6SEM_2D_QT_APPROXIMATE_TOOLS_H
