/*******************************************************************************
 *
 *		File "zeros.c"
 *
 * Files containing routines for finding zeros of functions
 * 
 * _ double Zbisection (double (*f)(double), double *v, double accuracy):
 * 		uses bisection method with final uncertainty specified by "accuracy".
 * 
 * _ double Zsecant (double (*f)(double), double *v, double accuracy):
 * 		uses secant method with final uncertainty specified by "accuracy".
 *
 ******************************************************************************/


#define ZEROS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Zbisection (double (*f)(double), double *v, double accuracy)
{
	if(f(v[0])*f(v[1])>=0)
	{
		printf("\nIt is impossible to use the bisection method in finding zeros of the function!\n\n");
		exit(EXIT_FAILURE);
	}
	
	double x;
	
	do{
		x = (v[0] + v[1])/2.0;
		if((f(v[0])*f(x))<0)
			v[1]=x;
		else if(f(x)==0)
			break;
		else
			v[0]=x;
	}while(fabs(v[0]-v[1]) > accuracy);
	
	return x;
}

double Zsecant (double (*f)(double), double *v, double accuracy)
{
	double x, f0, f1, fx;
	f0 = f(v[0]);
	f1 = f(v[1]);
	
	if(f0==f1 || (f0*f1>0))
	{
		printf("\nIt is impossible to use the secant method in finding zeros of the function!\n\n");
		exit(EXIT_FAILURE);
	}
	
	do{
		x = v[0] + (v[0]-v[1])*f0/(f1-f0);
		fx = f(x);
		
		if((f0*fx)<0)
		{	v[1]=x;   f1 = fx;   }
		
		else if(fx==0)
			break;
		else
			v[0]=x;
	}while(fabs(v[0]-v[1]) > accuracy);
	
	return x;
}
