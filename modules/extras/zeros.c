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
#include "extras.h"


double Zbisection (double (*f)(double), double *v, double accuracy)
{
	int i = 0;
	double x, f0, f1, fx;
	f0 = f(v[0]);
	f1 = f(v[1]);
	printf("\t\t x_min\t\t x_max\t\t f(x_min)\t f(x_max)\n");
	printf("Bisection %d :\t%e\t%e\t%e\t%e\n", i, v[0], v[1], f0, f1);
	
	if(f0==f1 || (f0*f1>0))
	{
		printf("\nIt is impossible to use the bisection method in finding the zeros of the function!\n\n");
		exit(EXIT_FAILURE);
	}
	
	do{
		i++;
		x = (v[0]+v[1])/2.0;
		fx = f(x);
		
		if((f0*fx)<0)
		{	v[1]=x;   f1 = fx;  }
		
		else if((f0*fx)>0)
		{	v[0]=x;   f0 = fx;  }

		else
			break;
			
		if(i==200){
			printf("\nBisection does not reach convergence!\n\n");
			exit(EXIT_FAILURE);
		}

		printf("Bisection %d :\t%e\t%e\t%e\t%e\n", i, v[0], v[1], f0, f1);

	}while(fabs(v[0]-v[1]) > accuracy);
	
	return x;
}

double Zsecant (double (*f)(double), double *v, double accuracy)
{
	int i=0;
	double x, f0, f1, fx;
	f0 = f(v[0]);
	f1 = f(v[1]);
	printf("\t\t x_min\t\t x_max\t\t f(x_min)\t f(x_max)\n");
	printf("Secant    %d :\t%e\t%e\t%e\t%e\n", i, v[0], v[1], f0, f1);
	if((f0==f1) || (f0*f1>0))
	{
		printf("\nIt is impossible to use the secant method in finding the zeros of the function!\n\n");
		exit(EXIT_FAILURE);
	}
	
	do{
		i++;
		x = v[0] + (v[0]-v[1])*f0/(f1-f0);
		fx = f(x);
		
		if((f0*fx)<0)
		{	v[1]=x;   f1 = fx;  }
		
		else if((f0*fx)>0)
		{	v[0]=x;   f0 = fx;  }

		else
			break;
		
		if(i==200){
			printf("\nSecant does not reach convergence!\n\n");
			exit(EXIT_FAILURE);
		}

		printf("Secant    %d :\t%e\t%e\t%e\t%e\n", i, v[0], v[1], f0, f1);

	}while(fabs(v[0]-v[1]) > accuracy);
	
	return x;
}
