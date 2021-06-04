#include"functions.h"

double eps = 0.00001;

double f0(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 1; }
double f1(double x, double y) { Q_UNUSED(y); return x; }
double f2(double x, double y) { Q_UNUSED(x); return y; }
double f3(double x, double y) { return x + y; }
double f4(double x, double y) { return sqrt(x*x + y*y); }
double f5(double x, double y) { return x*x + y*y; }
double f6(double x, double y) { return exp(x*x - y*y); }
double f7(double x, double y) { return 1/(25*(x*x + y*y) + 1); }

double Dxf0(double x, double y){ Q_UNUSED(x); Q_UNUSED(y); return 0; }
double Dxf1(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 1; }
double Dxf2(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double Dxf3(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 1; }
double Dxf4(double x, double y) 
{
	if (x*x + y*y > eps)
		return x/sqrt(x*x + y*y); 
	else return 0;
}
double Dxf5(double x, double y) { Q_UNUSED(y); return 2*x; }
double Dxf6(double x, double y) { return 2*x*exp(x*x - y*y); }
double Dxf7(double x, double y) { return -50*x/(25*(x*x + y*y) + 1)/(25*(x*x + y*y) + 1); }

double Dyf0(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double Dyf1(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double Dyf2(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 1; }
double Dyf3(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 1; }
double Dyf4(double x, double y) 
{
	if (x*x + y*y > eps)
		return y/sqrt(x*x + y*y); 
	else return 0;
}
double Dyf5(double x, double y) { Q_UNUSED(x); return 2*y; }
double Dyf6(double x, double y) { return -2*y*exp(x*x - y*y); }
double Dyf7(double x, double y) { return -50*y/(25*(x*x + y*y) + 1)/(25*(x*x + y*y) + 1); }

double DxDyf0(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double DxDyf1(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double DxDyf2(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double DxDyf3(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double DxDyf4(double x, double y) 
{
	if (x*x + y*y > eps)
		return -x*y/sqrt(x*x + y*y)/(x*x + y*y); 
	else return 0;
}
double DxDyf5(double x, double y) { Q_UNUSED(x); Q_UNUSED(y); return 0; }
double DxDyf6(double x, double y) { return -4*x*y*exp(x*x - y*y); }
double DxDyf7(double x, double y) { return 5000*y*x/(25*(x*x + y*y) + 1)/(25*(x*x + y*y) + 1)/(25*(x*x + y*y) + 1); }

double Pf(int i, int j, int NX, double x, double y, double* X, double* Y, double* coeffs)
{
	double res = 0;

	for (int k = 0; k < 4; ++k)
		for (int l = 0; l < 4; ++l)
			res += coeffs[(j*NX + i)*16 + 4*k + l] * W(i, j, k, l, x, y, X, Y);
			// res += G(i, j, k, l, NX, coeffs, Ax, Ay) * W(i, j, k, l, x, y, X, Y);

	return res;
}

double W(int i, int j, int k, int l, double x, double y, double* X, double* Y)
{
	return pow(x - X[i], k) * pow(y - Y[j], l);
}

void calculateCoeffs1(int NX, int NY, double* Ax, double* Ay, double* precoeffs, double* coeffs)
{
	double temp;

	for (int i = 0; i < NX; ++i)
		for (int j = 0; j < NY; ++j)
			for (int k = 0; k < 4; ++k)
				for (int l = 0; l < 4; ++l)
				{
					coeffs[(j*NX + i)*16 + k*4+l] = 0;
					for (int r = 0; r < 4; ++r)
					{
						temp = 0;
						for (int s = 0; s < 4; ++s)
							temp += precoeffs[(j*NX + i)*16 + r*4+s] * Ay[l*4+s];

						coeffs[(j*NX + i)*16 + k*4+l] += Ax[k*4+r] * temp;
					}
				}
}

// double L(double x1, double y1, double x2, double y2, double x, double y)
// {
// 	return (x - x1)*(y2 - y1) - (y - y1)*(x2 - x1);
// }

// double nf1(int ix, int iy, double *X, double *Y, double x, double y)
// {return L(X[ix*3+3], Y[iy*3], X[ix*3], Y[iy*3+3], x, y)/L(X[ix*3+3], Y[iy*3], X[ix*3], Y[iy*3+3], X[ix*3], Y[iy*3]);}
// double nf2(int ix, int iy, double *X, double *Y, double x, double y)
// {return L(X[ix*3], Y[iy*3], X[ix*3], Y[iy*3+3], x, y)/L(X[ix*3], Y[iy*3], X[ix*3], Y[iy*3+3], X[ix*3+3], Y[iy*3]);}
// double nf3(int ix, int iy, double *X, double *Y, double x, double y)
// {return L(X[ix*3], Y[iy*3], X[ix*3+3], Y[iy*3], x, y)/L(X[ix*3], Y[iy*3], X[ix*3+3], Y[iy*3], X[ix*3], Y[iy*3+3]);}

// double Pf1(int ix, int iy, double *X, double *Y, double *F, double x, double y)
// {
// 	return
// 		F[iy*3*((NX - 1)*3 + 1) + ix*3] * () + //1
// 		F[iy*3*((NX - 1)*3 + 1) + ix*3 + 3] * () + //2
// 		F[(iy*3 + 3)*((NX - 1)*3 + 1) + ix*3] * () + //3
// 		F[iy*3*((NX - 1)*3 + 1) + ix*3 + 1] * () + //4
// 		F[iy*3*((NX - 1)*3 + 1) + ix*3 + 2] * () + //5
// 		F[(iy*3 + 1)*((NX - 1)*3 + 1) + ix*3 + 2] * () + //6
// 		F[(iy*3 + 2)*((NX - 1)*3 + 1) + ix*3 + 1] * () + //7
// 		F[(iy*3 + 1)*((NX - 1)*3 + 1) + ix*3] * () + //8
// 		F[(iy*3 + 2)*((NX - 1)*3 + 1) + ix*3] * () + //9
// 		F[(iy*3 + 1)*((NX - 1)*3 + 1) + ix*3 + 1] * () + //10
// 		;
// }

// double RR(double f1, double f2, double x1, double x2) { return (f2 - f1)/(x2 - x1); }

//___________________________________

// void calculateCoeffs1(int N, double *X, double *F, double *A, double *D)
// {
// 	for (int i = 0; i < N - 1; ++i)
// 	{
// 		A[i*4] = F[i];
// 		A[i*4 + 1] = D[i];
// 		A[i*4 + 2] = RR(D[i], RR(F[i], F[i+1], X[i], X[i+1]), X[i], X[i+1]);
// 		A[i*4 + 3] = (D[i] + D[i+1] - 2 * RR(F[i], F[i+1], X[i], X[i+1])) / ((X[i+1] - X[i])*(X[i+1] - X[i]));
// 	}
// }

// double Pf1(int i, double* A, double x1, double x2, double x)
// {
// 	return A[4*i] + A[4*i+1]*(x - x1) + A[4*i+2]*(x - x1)*(x - x1) + A[4*i+3]*(x - x1)*(x - x1)*(x - x2);
// }

// //___________________________________

// double findD2(double f12, double f23)
// {
// 	if(f23 > 0)
// 	{
// 		if (f12 > 0)
// 			return min(fabs(f12), fabs(f23));
// 		else 
// 			return 0;
// 	}
// 	else
// 	{
// 		if (f12 < 0)
// 			return -min(fabs(f12), fabs(f23));
// 		else 
// 			return 0;
// 	}
// }

// void calculateCoeffs2(int N, double *X, double *F, double *A, double *D)
// {
// 	double f12, f23;

// 	D[0] = 0;
// 	D[N-1] = 0;

// 	for (int i = 1; i < N - 1; ++i)
// 	{
// 		f12 = RR(F[i-1], F[i], X[i-1], X[i]);
// 		f23 = RR(F[i], F[i+1], X[i], X[i+1]);

// 		D[i] = findD2(f12, f23);
// 	}

// 	for (int i = 1; i < N - 4; ++i)
// 	{
// 		A[i*4] = F[i+1];
// 		A[i*4 + 1] = D[i+1];
// 		A[i*4 + 2] = RR(D[i+1], RR(F[i+1], F[i+2], X[i+1], X[i+2]), X[i+1], X[i+2]);
// 		A[i*4 + 3] = (D[i+1] + D[i+2] - 2 * RR(F[i+1], F[i+2], X[i+1], X[i+2])) / ((X[i+2] - X[i+1])*(X[i+2] - X[i+1]));
// 	}

// 	A[0] = F[0];
// 	A[1] = RR(F[0], F[1], X[0], X[1]);
// 	A[2] = (F[2] - F[0] - A[1]*(X[2] - X[0])) / ((X[2] - X[0])*(X[2] - X[1]));
// 	A[3] = (D[2] - A[1] - A[2]*(2*X[2] - X[0] - X[1])) / ((X[2] - X[0])*(X[2] - X[1]));

// 	A[(N-4)*4] = F[N-1];
// 	A[(N-4)*4 + 1] = RR(F[N-1], F[N-2], X[N-1], X[N-2]);
// 	A[(N-4)*4 + 2] = (F[N-3] - F[N-1] - A[(N-4)*4 + 1]*(X[N-3] - X[N-1])) / ((X[N-3] - X[N-1])*(X[N-3] - X[N-2]));
// 	A[(N-4)*4 + 3] = (D[N-3] - A[(N-4)*4 + 1] - A[(N-4)*4 + 2]*(2*X[N-3] - X[N-1] - X[N-2])) / ((X[N-3] - X[N-1])*(X[N-3] - X[N-2]));

// }

// double Pf2(int i, double* A, double x1, double x2, double x3, double x)
// {
// 	return A[4*i] + A[4*i+1]*(x - x1) + A[4*i+2]*(x - x1)*(x - x2) + A[4*i+3]*(x - x1)*(x - x2)*(x - x3);
// }