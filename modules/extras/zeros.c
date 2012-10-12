/*******************************************************************************
 *
 *		File "zeros.c"
 *
 * Files containing routines for finding zeros of functions
 * 
 * _ double Zbisection (double (*f)(double), double *v, double precision):
 * 		uses bisection method with final uncertainty specified by "precision".
 * 
 * _ double Zsecant (double (*f)(double), double *v, double precision):
 * 		uses secant method with final uncertainty specified by "precision".
 *
 ******************************************************************************/


#define ZEROS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Zbisection (double (*f)(double), double *v, double precision)
{
	if(f(v[0])*f(v[1])>=0)
	{
		printf("\nIt is impossible to use the bisection method in finding zeros of the function!\n\n");
		exit(EXIT_SUCCESS);
	}
	
	double x = (v[0] + v[1])/2.0;
	if(fabs(v[0]-v[1])<precision)
		return x;
	
	else if(f(v[0])*f(x)<0)
		v[1] = x;
	else v[0] = x;
	
	Zbisection(f,v,precision);
}

double Zsecant (double (*f)(double), double *v, double precision)
{		
	if(f(v[0])==f(v[1]) || (f(v[0])*f(v[1])>0))
	{
		printf("\nIt is impossible to use the secant method in finding zeros of the function!\n\n");
		exit(EXIT_SUCCESS);
	}
	
	double x = v[0] + (v[0]-v[1])*f(v[0])/(f(v[1])-f(v[0]));
	if(fabs(v[0]-v[1])<precision)
		return x;
	
	else if(f(v[0])*f(x)<0)
		v[1] = x;
	else v[0] = x;
	
	Zsecant(f,v,precision);
}
