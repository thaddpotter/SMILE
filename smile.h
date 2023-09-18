#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Constants
const double pi = 3.1415926535;
const double dtor = pi/180;
const float sday = 86400;                 //Seconds per day
const float t0 = 3.2;                     //Pulse Midpoint
const double bigg = 6.667e-11;      //Grav Const
const double c = 2.99792458e8;      //Speed of light
const double ms = 1.988e30;         //Solar Mass

//Function Primitives
void shift_time(double *x1_0,double *y1_0, double *vx1_0,double *vy1_0,
            double *x2_0,double *y2_0, double *vx2_0,double *vy2_0,
            float ti, float m1, float m2, double rs);
void orbit(double x1_0, double y1_0, double vx1_0,double vy1_0,
            double x2_0, double y2_0, double vx2_0,double vy2_0,
            double x1[2048], double y1[2048], double vx1[2048],double vy1[2048],
            double x2[2048], double y2[2048], double vx2[2048],double vy2[2048], double dt, float m1, float m2, double rs);

double integrand(double u, double r, double t, double p, float gamma);
int big_integral(double pstar[2048], double u[2048], double init_curve[2048],float gamma);
int quadquad(double pstar[2048], double u[2048], double init_curve[2048], float gamma);
int piecequad(double pstar[2048], double u[2048], double init_curve[2048], float gamma);

int gen_light_curve(
            double times[2048], double light_curve[2048], float m1, float m2,
            float rstar, float ecc, float phi, float omega, float alpha,
            float gam, float g, float Periods, float ti);
