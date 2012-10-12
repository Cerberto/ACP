
#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#ifndef ODE2_C
extern void ev_K (double x1, double x2, double t, double dt, double (*f1)(double,double,double), \
	double (*f2)(double,double,double), double *K1, double *K2);
extern double deltaX (double *K);
#endif

#endif

