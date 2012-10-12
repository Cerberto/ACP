
/*******************************************************************************
 *
 *		File "methods.c"
 *
 * Metodi deterministici per il calcolo di integrali definiti.
 * Le routines sono:
 * _ ftrapezoid, dtrapezoid -> integrazione con metodo dei trapezi con numeri
 *		float e double rispettivamente;
 * _ fsimpson, dsimpson -> integrazione con metodo delle parabole con numeri
 *		float e double rispettivamente;
 * _ gaussianquad_5 -> integrazione con metodo delle quadrature gaussiane con
 *		polinomio di Legendre di gradi 5.
 *
 ******************************************************************************/

#define METHODS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integration.h"


float ftrapezoid (float (*f)(float), float a, float b)
{
	return (f(a) + f(b))*(b-a)/2;
}

double trapezoid (double (*f)(double), double a, double b)
{
	return (f(a) + f(b))*(b-a)/2;
}


float fsimpson (float (*f)(float), float a, float b)
{
	return (f(a) + 4*f((a+b)/2) + f(b))*(b-a)/6;
}


double simpson (double (*f)(double), double a, double b)
{
	return (f(a) + 4*f((a+b)/2) + f(b))*(b-a)/6;
}


double gaussianquad_5 (double (*f)(double), double a, double b)
{
	if(a == b)
		return 0;
	
	double S = 0;
	int i = 0;
	double roots[5];
	double weights[5];
	double x;
	
	/* Radici del polinomio di Legendre di gradi 5 nell'intervallo [-1,1] */
	roots[0] = 0;
	roots[1] = -0.53846931010568309;
	roots[2] = 0.53846931010568309;
	roots[3] = -0.90617984593866399;
	roots[4] = 0.90617984593866399;
	
	/* Pesi */
	weights[0] = 128.0/225.0;
	weights[1] = (322.0 +13.0*sqrt(70))/900.0;
	weights[2] = (322.0 +13.0*sqrt(70))/900.0;
	weights[3] = (322.0 -13.0*sqrt(70))/900.0;
	weights[4] = (322.0 -13.0*sqrt(70))/900.0;
	
	/* Calcolo dell'integrale */
	for(i = 0; i < 5; i++)
	{
		x = (b-a)/2 * roots[i] + (a+b)/2;
		S += weights[i] * f(x);
	}
	return (b-a)/2*S;
}

