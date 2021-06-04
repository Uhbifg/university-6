int sign(double x);
int method_init_1(int n, double *x, double *f_vals, double *additional_space, double *der);
double dd(double *f_vals, double *x, int i, int j);
int method_init_2(int n, double *x, double *f_vals, double *additional_space, double *add_nodes);
double w_j(double *f_vals, double *x, int j);
double max(double x1, double x2);


void matrix_product(double *a, double *b, int n, int m, int p, double* ans);
void get_a(double* ans, double h);
void method_1_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F);
double method_compute(double x, double y, double *x_vals, double *y_vals, int nx, int ny, double* F);
void transpose(double *a, int n);
void method_2_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F, double ax, double bx, double ay, double by,
                   double max_f, double x2, double y2, int p);

