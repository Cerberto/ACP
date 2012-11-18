
#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include "cluster.h"

#ifndef HARMOSC_C
extern double HOmetropolis (double *state, int state_dim);
extern double HOeAction (double *state, int state_dim);
/* Pay attention: replace this name to "autocorrelation" in harmosc_metropolis.c */
extern void HOautocorrelation (double *x, int t, double *v, int x_dim, int steps);
extern cluster DeltaE (cluster *A, cluster *B, cluster *C);
extern cluster MatrixElementX (cluster *DE, cluster *Corr, int t, int N);
extern cluster sqrt_jk (cluster *X);
#endif

#ifndef METROPOLIS_C
extern double autocorrelation (double *x, int t, int dim);
extern void cold_init (double *v, int dim);
extern void hot_init (double *v, int dim);
extern void metropolis (double (*P)(double *), double *state, int state_dim, double delta);
#endif

#endif
