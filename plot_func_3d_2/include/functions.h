#pragma once

#include<stdlib.h>
#include<math.h>
#include<string>

#include<QtGlobal>

using namespace std;

double f0(double, double);
double f1(double, double);
double f2(double, double);
double f3(double, double);
double f4(double, double);
double f5(double, double);
double f6(double, double);
double f7(double, double);

double Dxf0(double, double);
double Dxf1(double, double);
double Dxf2(double, double);
double Dxf3(double, double);
double Dxf4(double, double);
double Dxf5(double, double);
double Dxf6(double, double);
double Dxf7(double, double);

double Dyf0(double, double);
double Dyf1(double, double);
double Dyf2(double, double);
double Dyf3(double, double);
double Dyf4(double, double);
double Dyf5(double, double);
double Dyf6(double, double);
double Dyf7(double, double);

double DxDyf0(double, double);
double DxDyf1(double, double);
double DxDyf2(double, double);
double DxDyf3(double, double);
double DxDyf4(double, double);
double DxDyf5(double, double);
double DxDyf6(double, double);
double DxDyf7(double, double);

double Pf(int, int, int, double, double, double*, double*, double*);

double W(int, int, int, int, double, double, double*, double*);

void calculateCoeffs1(int, int, double*, double*, double*, double*);

// double L(double, double, double, double, double, double);

// double nf1(int, int, double*, double*, double, double);
// double nf2(int, int, double*, double*, double, double);
// double nf3(int, int, double*, double*, double, double);
// double nf4(int, int, double*, double*, double, double);
// double nf5(int, int, double*, double*, double, double);
// double nf6(int, int, double*, double*, double, double);
// double nf7(int, int, double*, double*, double, double);
// double nf8(int, int, double*, double*, double, double);
// double nf9(int, int, double*, double*, double, double);

// double Pf1(int, int, double*, double*, double*, double, double);

// double RR(double, double, double, double);

// void calculateCoeffs1(int, double*, double*, double*, double*);
// double Pf1(int, double*, double, double, double);


// double findD2(double, double);
// void calculateCoeffs2(int, double*, double*, double*, double*);
// double Pf2(int, double*, double, double, double, double);
