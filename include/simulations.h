
#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include "cluster.h"

#ifndef HARMOSC_C
extern double HOmetropolis(double *state, int state_dim);
extern double HOeAction(double *state, int state_dim);
extern void autocorrelation(double *x, int t, double *v, int x_dim, int steps);
extern void cold_init(double *v, int dim);
extern void hot_init(double *v, int dim);
extern cluster DeltaE(cluster *A, cluster *B, cluster *C);
extern cluster MatrixElementX(cluster *DE, cluster *Corr, int t, int N);
extern cluster sqrt_jk (cluster *X);
#endif

#endif
