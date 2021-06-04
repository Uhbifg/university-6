int sign(double x);
int method_init_1(int n, double *x, double *f_vals, double *a, double *additional_space, double *der);
double dd(double *f_vals, double *x, int i, int j);
int method_init_2(int n, double *x, double *f_vals, double *a, double *additional_space, double *add_nodes);
double w_j(double *f_vals, double *x, int j);
double max(double x1, double x2);


void matrix_product(double *a, double *b, int n, int m, int p, double* ans);
double* get_a(double h);
void method_1_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, auto F);
double method_compute(double x, double y, double *x_vals, double *y_vals, int nx, int ny, auto F);



