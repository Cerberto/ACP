#ifndef SCHRNUMEROV_H
#define SCHRNUMEROV_H

#define M 1.0

#ifndef NUMEROV_C
extern void evol (double (*V)(double), double (*S)(double), double r, double h, double * y, double E, int L);
extern void evol_back (double (*V)(double), double (*S)(double), double r, double h, double * y, double E, int L);
#endif


#endif
