int sign(double x);
int method_init_1(int n, double *x, double *f_vals, double *additional_space, double *der);
double dd(double *f_vals, double *x, int i, int j);
int method_init_2(int n, double *x, double *f_vals, double *additional_space, double *add_nodes);
double w_j(double *f_vals, double *x, int j);
double max(double x1, double x2);


void get_a(double* ans, double h);
void method_1_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F);
double method_compute(double x, double y, double *x_vals, double *y_vals, int nx, int ny, double* F);
void method_2_prod(double *fvals, int nx, int ny, double *x_vals, double *y_vals, int func_id, double* F, double ax, double bx, double ay, double by,
                   double max_f, double x2, double y2, int p);
int create_gamma(int nx, int ny, double*F, double *Fx, double *Fy, double *Fxy, double *fvals, double *y_vals, double *x_vals, int func_id);
double method_compute(double x, double y, double *x_vals, double *y_vals, int nx, int ny, double*  F, int i, int j);
double pow_d(double x, int p);
