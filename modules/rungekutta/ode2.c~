
/*******************************************************************************
 *
 *		File "ode2.c"
 *
 * Routines per la soluzione numerica di ODE lineari del secondo ordine con
 * condizioni iniziali con il metodo di Runge-Kutta al quarto ordine.
 *
 * Le funzioni sono:
 * _ ev_K -> assegnamento delle variabili ausiliarie K;
 * _ deltaX -> evoluzione della coordinata per un passo reticolare (restituisce
 *		lo "spostamento".
 *
 ******************************************************************************/

#define ODE2_C

#include <stdlib.h>
#include <math.h>

void ev_K (double x1, double x2, double t, double dt, double (*f1)(double,double,double), \
	double (*f2)(double,double,double), double *K1, double *K2)
{
	K1[0] = dt*f1(x1,x2,t);
	K2[0] = dt*f2(x1,x2,t);
	K1[1] = dt*f1(x1+0.5*K1[0], x2+0.5*K2[0], t+0.5*dt);
	K2[1] = dt*f2(x1+0.5*K1[0], x2+0.5*K2[0], t+0.5*dt);
	K1[2] = dt*f1(x1+0.5*K1[1], x2+0.5*K2[1], t+0.5*dt);
	K2[2] = dt*f2(x1+0.5*K1[1], x2+0.5*K2[1], t+0.5*dt);
	K1[3] = dt*f1(x1+K1[2], x2+K2[2], t+dt);
	K2[3] = dt*f2(x1+K1[2], x2+K2[2], t+dt);
}


double deltaX (double *K)
{
	return 1.0/6.0*(K[0] + K[3]) + 1.0/3.0*(K[1] + K[2]);
}
