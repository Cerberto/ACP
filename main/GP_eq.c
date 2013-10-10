
/**
 * 		File "GP_eq.c"
 * 
 * Self-consistent solution of the Gross-Pitaevskii equation
 * with the Numerov method in the case of central potential.
 * 
 * The sub-routines are:
 * _ init -> initialize an array with some given function;
 * _ potential -> sets the local potential present in the GP equation;
 * _ yRmax -> gives the value at the end of the mesh of the solution of
 * 		the Schrodinger-like equation calculated with the Numerov
 * 		method;
 * _ save / saveandprint -> save the best solution found with Numerov
 * 		within an array;
 * _ normalization -> compute the integral of a function (supposing
 * 		spherical symmetry) whose values are stored in an array.
 * 
 **/

#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "numerov.h"
#include "extras.h"

/**
 *		Constant parameters to set properly
 **/

/* "Normalization" constant */
#define N 1.0e-2

/* length and energy steps */
#define H 5.0e-4
#define Dmu 0.01

/* accuracy for bisection and secant methods */
#define B_ACC 1.0e-3
#define S_ACC 1.0e-6

/* Mixing constant (for the new density) */
#define BETA 0.983

/* Maximum radius */
#define RMAX 6.95

/* Width of initial density */
#define SIGMA 1.5

double PI = 4*atan(1.0);


/**
 * 		External variables and sub-routines' definition
 **/

/* Number of particles and scattering length */
int NPAR;
double ALPHA;

/* density and number of points of the mesh */
double *rho; int DIM;

/* array for the solution */
double *y;

void init (double *v, int dim);
double potential (double x);
double yRmax (double E);		/* Function that computes the value of the wave function at the end of the mesh */

/* functions for saving (and printing) solutions */
void save (double E, double *x, int dim);
void saveandprint (double E, double *x, int dim, char *filename);


double normalization(double *v, int dim);

/**
 * 		Main routine
 **/
int main (int argc, char *argv[])
{
	if(argc<4)
	{
		printf("\nMissing input parameters: please specify precision needed for the eigenvalue, number of particles and scattering length.\n\n");
		exit(EXIT_FAILURE);
	}
	
	int counter = 0;
	
	/* Chemical potential mu and its range of variation */
	double mu;
	double muMIN, epsilon, EV, muOLD;
	epsilon	= atof(argv[1]);	// accuracy of chemical potential
	
	NPAR	= atoi(argv[2]);
	ALPHA	= atof(argv[3]);
	/* set the starting value of the chemical potential */
		if(ALPHA < 0)
			muMIN = 4*PI*ALPHA + 1.0e-05;
		else
			muMIN = 1.0e-05;
	EV		= -100.0;
	muOLD	= -1000.0;
	double temp		= 0;
	double mutemp	= 0;
	
	DIM = (int)(RMAX/H);
	rho = malloc(DIM*sizeof(double));
	y = malloc(DIM*sizeof(double));
	
	/* Extremes of the range where searching for the zero */
	double *ye;
	ye = malloc(2*sizeof(double));
	
	/* Preparing files and folder to save outputs */
	char *path, *out_file, *command;
		path		= malloc(100*sizeof(char));
		out_file	= malloc(100*sizeof(char));
		command		= malloc(100*sizeof(char));
	FILE *output_ev, *mu_vs_alpha;
	sprintf(path, "GP/%d", NPAR);
	sprintf(out_file, "%s/mu_vs_alpha.dat", path);
		mu_vs_alpha = fopen(out_file,"a");
	sprintf(command, "mkdir %s/%.4lf", path, ALPHA);
		system(command);
	sprintf(out_file, "/%.4lf/eigenvalues.dat", ALPHA);
	strcat(path, out_file);
		output_ev	= fopen(path, "w");
	
	printf("\nScattering length = %lf.\n", ALPHA);
	printf("\n%d particles.\n\n", NPAR);
	
	init(rho, DIM);

	/* Change chemical potential until convergence is reached 
	 * and GP equation is considered to be solved
	 */
	do{
		counter++;
		mu = muMIN;
		
		/*
		 * For each value of the chemical potential we find the solution
		 * of the differential equation
		 */
		while(1)
		{
			/* solution of the equation for a given value of mu */
			yRmax(mu);
		
			/* compare value of solution at the end of the mesh with the previous result
			 * and if it has opposite sign it looks for zero between the last two
			 * values of mu */
			if((temp*y[DIM-1])<0)
			{
				ye[0]=mutemp;
				ye[1]=mu;
				printf("Before search:\t%e\t%e\t%e\t%e\n\n", mutemp, mu, temp, y[DIM-1]);
				muOLD = EV;
				Zbisection(yRmax, ye, B_ACC);		// bisection method up to a certain precision
				EV = Zsecant(yRmax, ye, S_ACC);		// refine search with secant method
				save(EV, y, DIM);
				
				printf("\n%d ... OLD mu = %.11e, \t NEW mu = %.11e\n\n", counter, muOLD, EV);
				fprintf(output_ev, "%d\t%.11e\n", counter, EV);		// print values of the chemical potential as function of iterations
				printf("scattering length = %.4lf\n", ALPHA);
				printf("normalization     = %.2lf\n\n", normalization(rho,DIM));
				break;
			}
			temp = y[DIM-1];
			mutemp = mu;
			mu += Dmu;
		}
	}while(fabs(muOLD-EV) > epsilon);
	
	/* Set output files and save the final solution */
	sprintf(path, "GP/%d/%.4lf", NPAR, ALPHA);
		/* Print on terminal the final results */
		printf("\n\nFinal wave function and density printed on \"%s/ultimate.dat\"", path);
		printf("\nChemical potentials as function of iteration printed on \"%s/eigenvalues.dat\"", path);
		printf("\nFinal chemical potential = %.4lf\n\n", EV);
		fprintf(mu_vs_alpha, "%.4lf\t%.6e\n", ALPHA, EV);
	sprintf(out_file, "/ultimate.dat");
	strcat(path, out_file);
	saveandprint(EV, y, DIM, path);
	
	fclose(output_ev);
	fclose(mu_vs_alpha);
	
	exit(EXIT_SUCCESS);
}


/**
 * 		Subroutines
 **/

void init (double *v, int dim)
{
	int i;
	for(i=0; i<dim; i++)
//		v[i] = 1.0/((i*H+H)*(i*H+H) + SIGMA*SIGMA);
//		v[i] = 1.0/(exp(-20*(i*H*H-SIGMA)) + 1.0);
		v[i] = exp(-(i*H*H)*(i*H*H));
	double S = normalization(v,dim);
	for(i=0; i<dim; i++)
		v[i] *= (NPAR/S);
}


double potential (double x)
{
	return 0.5*x*x;
}


double yRmax (double E)
{
	y[0] = H*N;	y[1] = 2*H*N;
	
	evol_GP(potential, rho, ALPHA, y, DIM, H, E);
	
	return y[DIM-1];
}


void saveandprint (double E, double *x, int dim, char *filename)
{
	int i;
	double S=0;
	FILE *output;
		output = fopen(filename, "w");
	x[0] = N*H;	x[1] = 2*H*N;
	evol_GP(potential, rho, ALPHA, x, DIM, H, E);
	
	/* Normalization of wave function */
	
	/*****************************************************************
	 * calcolo integrale di y^2 = r^2 * R^2 da 0 a infinito:
	 * questo deve essere uguale al numero di particelle
	 *****************************************************************/
	for(i=dim-1; i>-1; i--)
		S += x[i]*x[i];
	S *= H;
	
	fprintf(output, "#\n#Columns are:\n");
	fprintf(output, "#  Radius\t\t Wavefunction\t\t Density\n#\n");	 
	for(i=0; i<DIM; i++)
	{
		rho[i] = x[i]/(i*H+H);	rho[i] *= (NPAR*rho[i]/(4*PI*S));
		fprintf(output,"%.11e\t%.11e\t%.11e\n", (i+1)*H, x[i], rho[i]);
	}
	
	fclose(output);
}


void save (double E, double *x, int dim)
{
	int i;
	double S=0;
	double new;
	x[0] = N*H;	x[1] = 2*H*N;
	evol_GP(potential, rho, ALPHA, x, DIM, H, E);
	
	/* Normalization of wave function */
	for(i=dim-1; i>-1; i--)
		S += x[i]*x[i];
	S *= H;
	
	for(i=0; i<dim; i++)
	{
		new = x[i]/(i*H+H);	new *= (NPAR*new/(4*PI*S));
		rho[i] = BETA*rho[i] + (1-BETA)*new;
	}
}

/* compute the normalization of the (SPHERICALLY SYMMETRIC) DENSITY,
 * so the 4*pi is there */
double normalization (double *v, int dim)
{
	int i;
	double S = 0;
	for(i=dim-1; i>-1; i--)
		S += v[i]*(i*H+H)*(i*H+H);
	return 4*PI*S*H;
}

