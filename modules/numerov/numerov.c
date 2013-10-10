/*******************************************************************************
 *
 *	File "numerov2.c"
 *
 ******************************************************************************/


#define NUMEROV_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerov.h"

static double PI = 4*atan(1.0);

static double Ksq (double x, double (*V)(double), double E, int L)
{
	if(L==0)
		return 2*M*(E - V(x));

	else if(L<0)
	{	printf("\nL<0; invalid parameter.\n\n");
		exit(EXIT_FAILURE);
	}

	else
		return 2*M*(E - L*(L+1)/(2*M*x*x) - V(x));
}


void evol (double (*V)(double), double (*S)(double), double r, double h, double *y, double E, int L)
{
	double Know, Kprev, Knext;
		Know	= Ksq(r,V,E,L);
		Kprev	= Ksq(r-h,V,E,L);
		Knext	= Ksq(r+h,V,E,L);
	
	 y[2] = (h*h/12.0*(S(r+h) + 10*S(r) + S(r-h)) + y[1]*(2 - 5.0/6.0*Know*h*h) - y[0]*(1 + h*h/12.0*Kprev))/(1 + h*h/12.0*Knext);
	 
	 y[0] = y[1];  y[1] = y[2];
}


void evol_back (double (*V)(double), double (*S)(double), double r, double h, double *y, double E, int L)
{
	double Know, Kprev, Knext;
		Know	= Ksq(r,V,E,L);
		Kprev	= Ksq(r-h,V,E,L);
		Knext	= Ksq(r+h,V,E,L);
	
	 y[0] = (h*h/12.0*(S(r+h) + 10*S(r) + S(r-h)) + y[1]*(2 - 5.0/6.0*Know*h*h) - y[0]*(1 + h*h/12.0*Knext))/(1 + h*h/12.0*Kprev);
	 
	 y[2] = y[1];  y[1] = y[0];
}


void evol_GP (double (*V)(double), double *rho, double alpha, double *y, int dim, double h, double mu)
{
	int i;
	double Know, Kprev, Knext;
	
	for(i=1; i<dim-1; i++)
	{
		Kprev	= 2*(mu - 4*PI*alpha*rho[i-1] - V(i*h-h));
		Know	= 2*(mu - 4*PI*alpha*rho[i] - V(i*h));
		Knext	= 2*(mu - 4*PI*alpha*rho[i+1] - V(i*h+h));
	
		y[i+1] = (y[i]*(2 - 5.0/6.0*Know*h*h) - y[i-1]*(1 + h*h/12.0*Kprev))/(1 + h*h/12.0*Knext);
	}
}


