
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
#include <math.h>
#include "schrnumerov.h"
#include "extras.h"


/* "Normalization" constant */
#define N 1.0e-2

/* Increment and number of steps */
//#define EMIN 3.5
#define DE 0.007
#define ESTEPS 1000

/* Step length */
#define H 1.0e-5

/* Maximum radius */
//#define RMAX 5.0
double RMAX = 5.0;

/* Proper frequency */
#define Omega 1.0

/* #define A 1.0
 * #define mu 1.0
 * #define R 1.0
 */

double *y;
int L;

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
double yRmax (double E)
{
	double r;
	int j;
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
	return y[1];
}


void print_solution (double E, double *x, char *filename)
{
	double r;
	int j;
	FILE *output;
		output = fopen(filename, "w");
	x[0] = N;	x[2] = 0;	x[1] = N;
	for(j=0; j<L+1; j++)
	{
		x[0]*=H;
		x[1]*=(2*H);
	}
	
	r = 2*H;
	while(r<RMAX)
	{
		evol(V, null, r, H, x, E, L);
		fprintf(output,"%.11e\t%.11e\n", r, x[1]);
		r += H;
	};
	
	fclose(output);
	
	/* METTERE QUI UNA CHIAMATA A SYSTEM PER STAMPARE SOLUZIONE CON GNUPLOT */
}


int main (int argc, char *argv[])
{
	if(argc<3)
	{
		printf("\nMissing input parameters: please specify L (angular momentum) and E.\n\n");
		exit(EXIT_FAILURE);
	}
	double EMIN;
		EMIN = atof(argv[2]);
	double E = EMIN;	/* energy */
	double EV;			/* energy eigenvalue */
	double temp = 0;
	double Etemp = 0;
	int i, zeros;
	L = atoi(argv[1]);
	
	/* array in which are going to be saved the 3 points used in the algorithm
	 * and their initialization */
	y = malloc(3*sizeof(double));
	
	double *ye, *X;
	ye = malloc(2*sizeof(double));
	X = malloc(3*sizeof(double));
	
	char *out_file;
		out_file = malloc(100*sizeof(char));
		
	
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
		RMAX = 2.3*sqrt(E + sqrt(E*E -L*L -L));
		yRmax(E);
		fprintf(outboundary, "%e\t%e\n", E, y[1]);
		
		printf("%.11e\t%.11e\t%.11e\n", E, temp, y[1]);
		fflush(stdout);
		
		if((temp*y[1])<0)
		{
			zeros++;
			ye[0]=Etemp;
			ye[1]=E;
			printf("%e\t%e\n", ye[0], ye[1]);
				fflush(stdout);
			Zbisection(yRmax, ye, 1.0e-4);
			EV = Zsecant(yRmax, ye, 1.0e-8);
			sprintf(out_file, "isotropic_HO/solution_%d_%.7lf.dat", L, EV);
			print_solution(EV, X, out_file);
		}
		temp = y[1];
		Etemp = E;
		E += DE;
	}
	
	fclose(outboundary);
	
	exit(EXIT_SUCCESS);
}
