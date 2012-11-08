

#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "numerov.h"

#include "extras.h"

/* "Normalization" constant */
#define N 1.0e-2

/* length and energy steps */
#define H 1.0e-5
#define Dmu 0.007

/* Interaction constant */
#define ALPHA 0.1

/* Maximum radius */
#define RMAX 5.0

double PI = 4*atan(1.0);

/* density and number of points of the mesh */
double *rho; int DIM;

/* array for the solution */
double *y;

void gauss_init (double *v, int dim);
double potential (double x);
double yRmax (double E);		/* Function that computes the value of the wave function at the end of the mesh */
/* functions for saving (and printing) solutions */
void save (double E, double *x, int dim);
void saveandprint (double E, double *x, int dim, char *filename);


int main (int argc, char *argv[])
{
	if(argc<2)
	{
		printf("\nMissing input parameters: please specify precision needed for the eigenvalue.\n\n");
		exit(EXIT_FAILURE);
	}
	
	int counter = 0;
	double mu;
	double muMIN, epsilon, EV, muOLD;
		if(ALPHA < 0)
			muMIN = 4*PI*ALPHA + 1.0e-05;
		else
			muMIN = 1.0e-05;
	epsilon	= atof(argv[1]);
	EV		= -100.0;
	muOLD	= -1000.0;
	double temp		= 0;
	double mutemp	= 0;
	
	DIM = (int)(RMAX/H);
	rho = malloc(DIM*sizeof(double));
	y = malloc(DIM*sizeof(double));
	
	double *ye;
	ye = malloc(2*sizeof(double));
	
	char *path, *out_file, *command;
		path		= malloc(100*sizeof(char));
		out_file	= malloc(100*sizeof(char));
		command		= malloc(100*sizeof(char));
		
	sprintf(path, "GP/GP_%.4lf", ALPHA);
	sprintf(command, "mkdir %s", path);
		system(command);
	sprintf(out_file, "/eigenvalues_%.4lf.dat", ALPHA);
	strcat(path, out_file);
		FILE *output_ev;
		output_ev = fopen(path, "w");
	
	gauss_init(rho, DIM);
	
	do{
		counter++;
		mu = muMIN;
		while(1)
		{
			/* solution of the equation for a given value of mu */
			yRmax(mu);
		
			/* compares value of solution at the end of the mesh with the previous result
			 * and if it has opposite sign it looks for zero between the last two
			 * values of mu */
			if((temp*y[DIM-1])<0)
			{
				ye[0]=mutemp;
				ye[1]=mu;
				muOLD = EV;
				Zbisection(yRmax, ye, 1.0e-4);		// bisection method up to a certain precision
				EV = Zsecant(yRmax, ye, 1.0e-8);	// refine search with secant method
				printf("%.9e\n", EV);
/*				sprintf(path, "GP/GP_%.4lf", ALPHA);
				sprintf(out_file, "/solution_%.4lf_%.7lf.dat", ALPHA, EV);
				strcat(path, out_file);
				saveandprint(EV, y, DIM, path);*/
				save(EV, y, DIM);
				printf("%d ... muOLD = %.11e, \t EV = %.11e\n", counter, muOLD, EV);
				fprintf(output_ev, "%d\t%.11e\n", counter, EV);		// print values of mu's as function of iterations
				break;
			}
			temp = y[DIM-1];
			mutemp = mu;
			mu += Dmu;
		}
	}while(fabs(muOLD-EV) > epsilon);
	
	sprintf(path, "GP/GP_%.4lf", ALPHA);
	sprintf(out_file, "/ultimate_%.4lf_%.7lf.dat", ALPHA, EV);
	strcat(path, out_file);
	saveandprint(EV, y, DIM, path);
	
	exit(EXIT_SUCCESS);
}


void gauss_init (double *v, int dim)
{
	int i;
	for(i=0; i<dim; i++)
		v[i] = exp(i*i*H*H);
}


double potential (double x)
{
	return 0.5*x*x;
}


double yRmax (double E)
{
	y[0] = H;	y[1] = 2*H;
	
	evol_GP(potential, rho, ALPHA, y, DIM, H, E);
	
	return y[DIM-1];
}


void saveandprint (double E, double *x, int dim, char *filename)
{
	int i;
	FILE *output;
		output = fopen(filename, "w");
	x[0] = H;	x[1] = 2*H;
	evol_GP(potential, rho, ALPHA, x, DIM, H, E);	
	for(i=1; i<DIM; i++)		/* ATTENZIONE AGLI INDICI!! FALLI PARTIRE GIUSTI!! */
	{
		rho[i] = x[i]/(double)(i*H);	rho[i] *= rho[i];
		fprintf(output,"%.11e\t%.11e\t%.11e\n", i*H, x[i], rho[i]);
	};
	
	fclose(output);
}

void save (double E, double *x, int dim)
{
	int i;
	x[0] = H;	x[1] = 2*H;
	evol_GP(potential, rho, ALPHA, x, DIM, H, E);	
	for(i=1; i<DIM; i++)
	{	rho[i] = x[i]/(double)(i*H);	rho[i] *= rho[i];	};
}

