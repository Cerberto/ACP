#ifndef NUMEROV_H
#define NUMEROV_H

#define M 1.0

#ifndef NUMEROV_C
extern void evol (double (*V)(double), double (*S)(double), double r, double h, double * y, double E, int L);
extern void evol_back (double (*V)(double), double (*S)(double), double r, double h, double * y, double E, int L);
void evol_GP (double (*V)(double), double *rho, double alpha, double *y, int dim, double h, double mu);
#endif


#endif
