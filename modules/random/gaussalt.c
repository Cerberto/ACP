/*******************************************************************************
*
*		File "gausswild.c"
* 
* Generation of Gaussian random numbers
*
* The externally accessible functions are
*
*   void wgauss_dble(double rd[],int n)
*     Generates n double-precision Gaussian random numbers x with distribution
*     proportional to exp(-x^2/(2*sigma^2)) and assigns them to rd[0],..,rd[n-1]
* 
* 	double wgaussdistr (double x, double sigma)
* 	  Function corresponding to the distribution according to which gauss and
* 	  gauss_dble generate random numbers, namely a gaussian distribution
* 	  with variance sigma^2.
*
*   void stdgauss_dble(double rd[],int n)
*     Generates n double-precision Gaussian random numbers x with distribution
*     proportional to exp(-x^2/2) and assigns them to rd[0],..,rd[n-1]
* 
* 	double wgaussdistr (double x, double sigma)
* 	  Function corresponding to the distribution according to which gauss and
* 	  gauss_dble generate random numbers, namely a gaussian distribution
* 	  with variance sigma^2.
*
*******************************************************************************/

#define GAUSSALT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"

#define PI 3.141592653589793



/*
 * Distribuzione gaussiana con sigma definita dall'esterno
 */

void wgauss_dble(double rd[],int n, double sigma)
{
	if (sigma <= 0)
	{
		printf("\nSigma non positiva!!\n\n");
		exit(EXIT_FAILURE);
	}
	
	else
	{
		int k;
		double u[2];
		double x1,x2,rho,y1,y2;

		for (k=0;k<n;)
		{
			ranlxd(u,2);
			x1=u[0];
			x2=u[1];

			rho=-2*sigma*sigma*log(1.0-x1);
			rho=sqrt(rho);
			x2*=2.0*PI;
			y1=rho*sin(x2);
			y2=rho*cos(x2);
      
			rd[k++]=y1;
			
			if (k<n)
				rd[k++]=y2;
		}
	}
}


double wgaussdistr(double x, double sigma)
{
	if (sigma <= 0)
	{
		printf("\nSigma non positiva!!\n\n");
		exit(EXIT_FAILURE);
	}
	
	else
		return exp(-x*x/2.0/sigma/sigma)/sqrt(2*PI*sigma*sigma);
}




/*
 * Distribuzione gaussiana standard
 */
 
void stdgauss_dble(double rd[],int n)
{
	int k;
	double u[2];
	double x1,x2,rho,y1,y2;

	for (k=0;k<n;)
	{
		ranlxd(u,2);
		x1=u[0];
		x2=u[1];

		rho=-2.0*log(1.0-x1);
		rho=sqrt(rho);
		x2*=2.0*PI;
		y1=rho*sin(x2);
		y2=rho*cos(x2);
      
		rd[k++]=y1;
			
		if (k<n)
			rd[k++]=y2;
	}
}


double stdgaussdistr(double x)
{
	return exp(-x*x/2.0)/sqrt(2.0*PI);
}

