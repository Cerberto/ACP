

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
#define ALPHA 0.01

/* Maximum radius */
#define RMAX 5.0

double PI = 4*atan(1.0);

/* density and numbe of points of the mesh */
double *rho; int DIM;

double *y;


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


void evol_GP (double (*V)(double), double *rho, double *y, int rho_dim, int i, double h, double mu)
{
	if(i<1)
	{
		printf("Error: passed index < 1 to \"evol_GP\" routine.\n\n");
		exit(EXIT_FAILURE);
	}
	
	double Know, Kprev, Knext;
		Kprev	= 2*(mu - 4*PI*ALPHA*rho[i-1] - V(i*h-h));
		Know	= 2*(mu - 4*PI*ALPHA*rho[i] - V(i*h));
		Knext	= 2*(mu - 4*PI*ALPHA*rho[i+1] - V(i*h+h));
	
	 y[2] = (y[1]*(2 - 5.0/6.0*Know*h*h) - y[0]*(1 + h*h/12.0*Kprev))/(1 + h*h/12.0*Knext);
	 
	 y[0] = y[1];  y[1] = y[2];
}


/* Function that computes the value of the wave function at the end of the mesh */
double yRmax (double E)
{
	int i;
	double r;
	y[0] = H;	y[2] = 0;	y[1] = 2*H;
	
	r = 2*H;
	for(i=1; i<DIM; i++)		/* ATTENZIONE AGLI INDICI!! FALLI PARTIRE GIUSTI!! */
	{
		evol_GP(potential, rho, y, DIM, i, H, E);
		r += H;
	};
	return y[1];
}


void saveandprint (double E, double *x, int dim, char *filename)
{
	double r;
	int i;
	FILE *output;
		output = fopen(filename, "w");
	x[0] = H;	x[2] = 0;	x[1] = 2*H;
	rho[0] = 1.0;
	r = 2*H;
	for(i=1; i<DIM; i++)		/* ATTENZIONE AGLI INDICI!! FALLI PARTIRE GIUSTI!! */
	{
		evol_GP(potential, rho, x, DIM, i, H, E);
		pippo();
		fprintf(output,"%.11e\t%.11e\n", r, x[1]);
		pippo();
		rho[i] = x[1]*x[1];		/* ATTENZIONE QUI!! COME E' RHO?!?! */
		r += H;
	};
	
	fclose(output);
}



int main (int argc, char *argv[])
{
	if(argc<3)
	{
		printf("\nMissing input parameters: please specify maximum value of mu and precision needed on its evaluation.\n\n");
		exit(EXIT_FAILURE);
	}
	
	double mu;
	double muMIN, muMAX, epsilon, EV, muOLD;
		if(ALPHA < 0)
			muMIN = 4*PI*ALPHA + 1.0e-05;
		else
			muMIN = 1.0e-05;
	muMAX	= atof(argv[1]);
	epsilon	= atof(argv[2]);
	EV		= 0;
	muOLD	= 0;
	double temp		= 0;
	double mutemp	= 0;
	
	DIM = (int)(RMAX/H);
	rho = malloc(DIM*sizeof(double));
	y = malloc(3*sizeof(double));
	
	double *ye, *X;
	ye = malloc(2*sizeof(double));
	X = malloc(3*sizeof(double));
	
	char *out_file;
		out_file = malloc(100*sizeof(char));

	gauss_init(rho, DIM);
	
	do{
		mu = muMIN;
		while(mu < muMAX)
		{
			/* evolve il sistema per un certo valore di mu */
			yRmax(mu);
		
			/* confronta col precedente e se discordi cerca lo 0 */
			/* trovato lo zero (EV) lo si salva in muOLD e si esce dal ciclo */
			if((temp*y[1])<0)
			{
				ye[0]=mutemp;
				ye[1]=mu;
				muOLD = EV;
				Zbisection(yRmax, ye, 1.0e-4);
				EV = Zsecant(yRmax, ye, 1.0e-8);
				printf("%.9e\n", EV);
//				sprintf(out_file, "GP/solution_%.7lf.dat", EV); 
				saveandprint(EV, X, DIM, out_file);
				break;
			}
			temp = y[1];
			mutemp = mu;
			mu += Dmu;
		}
	}while(fabs(muOLD-EV) > epsilon);
	
	sprintf(out_file, "GP/ultimate_%.2lf_%.5e.dat", ALPHA, EV);
	saveandprint(EV, X, DIM, out_file);
	
	exit(EXIT_SUCCESS);
}
