
/*******************************************************************************
 *
 *	Program for solving Schrodinger equation in the case of isotropic harmonic
 * oscillator with Numerov method.
 *	Routines are defined in ../modules/numerov/numerov2.c
 *
 ******************************************************************************/
 
#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "schrnumerov.h"
#include "extras.h"

/* Maximum radius */
#define RMAX 5.0

/* "Normalization" constant */
#define N 1.0e-2

/* Increment and number of steps */
#define EMIN 3.5
#define DE 0.01
#define ESTEPS 500

/* Step length */
#define H 1.0e-5

/* Proper frequency */
#define Omega 1.0

/* #define A 1.0
 * #define mu 1.0
 * #define R 1.0
 */

/* Potential */
double V (double r)
{	
	return 0.5*M*Omega*Omega*r*r;
	//return -A/(1+ exp(mu*(r-R)));
}


/* Null function */
double null (double x)
{
	return 0;
}


/* Function that computes the value of the wave function at the end of the mesh */
double yRmax (double *y, double E, int L)
{
	y[0] = N;	y[2] = 0;	y[1] = N;
	for(j=0; j<L+1; j++)
	{
		y[0]*=H;
		y[1]*=(2*H);
	}
	
	r = 2*H;
	while(r<RMAX)
	{
		evol(V, null, r, H, y, E, L);
		r += H;
	};
}


int main (int argc, char *argv[])
{
	if(argc<2)
	{
		printf("\nMissing input parameters: please specify L (angular momentum)\n\n");
		exit(EXIT_FAILURE);
	}
	
	double r;

	double E = EMIN;
	int i, j, L;
	int zeros;
	L = atoi(argv[1]);
	
	/* array in which are going to be saved the 3 points used in the algorithm
	 * and their initialization */
	double y[3];
	
	/* array in which will be saved two values of E and the correspondent values of the solution at the last point of the mesh;
	 * this will be used in finding the eigenvalues
	 **/
	double Ey[4];
	
	char *out_file;
		out_file = malloc(100*sizeof(char));
	/*FILE *output[ESTEPS];
	for(i=0; i<ESTEPS; i++)
	{
		sprintf(out_file, "isotropic_HO/solution_%.5lf.dat", EMIN+i*DE);
		//sprintf(out_file, "wood_saxon/solution_%.5lf.dat", EMIN+i*DE);
		output[i] = fopen(out_file, "w");
	}*/
	FILE *outboundary;
		sprintf(out_file, "isotropic_HO/yRmaxE_%d.dat", L);
		//sprintf(out_file, "wood_saxon/yRmaxE_%d.dat", L);
	outboundary	= fopen(out_file,"w");
	fprintf(outboundary,
		"#\n# Value of the radial part of the isotropic harmonic oscillator's solution\n# at the last point of the mesh as a function of the energy parameter.\n# Angular momentum (absolute value) quantum number L = %d\n#\n", L);
	/* fprintf(outboundary,
		"#\n# Value of the radial part of the Wood-Saxon well's solution\n# at the last point of the mesh as a function of the energy parameter.\n# Angular momentum (absolute value) quantum number L = %d\n#\n", L); */
	
	for(i=0; i<ESTEPS; i++)
	{
		printf("\tStep %d of %d\n", i+1, ESTEPS);
		fflush(stdout);
		yRmax(y, E, L);
		fprintf(outboundary, "%e\t%e\n", E, y[1]);
		if(temp*y[1]<0)
		{
			zeros++;
			Ey[0]=Etemp;
			Ey[1]=E;
			Ey[2]=ytemp;
			Ey[3]=y[1];
			eigenvalue( /*...*/);
		}
		temp = y[1];
		r = H;
		Etemp = E;
		E += DE;
	}
	
	fclose(outboundary);
/*	for(i=0; i<ESTEPS; i++)
		fclose(output[i]);*/
	
	exit(EXIT_SUCCESS);
}
