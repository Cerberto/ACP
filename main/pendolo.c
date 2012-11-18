
/*******************************************************************************
 *
 *		File "pendolo.c"
 *
 *	Soluzione numerica dell'equazione del pendolo caotico
 * con il metodo di Runge-Kutta al IV ordine.
 *
 ******************************************************************************/
 
#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "extras.h"
#include "rungekutta.h"

/* Tempo di evoluzione del sistema */
#define T 200.0

/* Coefficiente di smorzamento */
#define Q 0.2

/* Ampiezza e frequenza della forzante esterna */
#define B 0.05
#define WF 0.1


/* funzione velocità */
double func1 (double x1, double x2, double  t)
{
	return x2;
}


/* funzione accelerazione */
double func2 (double x1, double x2, double t)
{
	return (- sin(x1) - Q*x2 + B*cos(WF*t));
}


int main (int argc, char *argv[])
{
	if(argc < 3)
	{
		printf("\nNumero di parametri insufficiente!\n");
		printf("Inserire:\n");
		printf("\t 1_ posizione iniziale;\n");
		printf("\t 2_ velocità iniziale;\n");
		printf("\t 3_ passo temporale -default: 0.01- (opzionale);\n");
		printf("\t 4_ file output -default:\"dati.dat\"- (opzionale).\n\n");
		exit(EXIT_FAILURE);
	}
	
	int step, STEPS;
	double H = 0.01;
	double t;
	double u;
	double X1, X2, Z1, Z2;
	double *K1, *K2, *J1, *J2;
	double Err1, Err2;
	double temp1[2], temp2[2];
	double (*velocity)(double, double, double);
	double (*acceleration)(double, double, double);

	char *file, *path;
	file = malloc(100*sizeof(char));
	path = malloc(100*sizeof(char));
		sprintf(path, "rungekutta/");
		sprintf(file, "dati.dat");
	STEPS = T/H;
	if((STEPS%2)!=0)
		printf("\nImpossibile calcolare l'errore\n\n");
	
	velocity		= func1;
	acceleration	= func2;
	
	// Condizioni iniziali
	X1 = atof(argv[1]);		// posizione
	X2 = atof(argv[2]);		// velocità
	Z1 = X1;
	Z2 = X2;
	
	if(argc == 4)
		H = atof(argv[3]);
	if(argc >= 5)
		file = argv[4];
	
	strcat(path, file);
	
	FILE *output, *errors;
	output = fopen(path,"w");
		fprintf(output, "#\n# Dati dell'evoluzione temporale del pendolo caotico;\n");
		fprintf(output, "# le colonne sono:\n");
		fprintf(output,	"# 1_tempo,  2_angolo (rispetto alla verticale),  3_velocità angolare.\n#\n");
	errors = fopen("rungekutta/errors.dat", "a");
/*		fprintf(errors, "#\n# Errore calcolato come differenza tra evoluti x_{n+2} con passo H e 2H\n");
		fprintf(errors, "# al variare di H; media su tutti i passi\n");
		fprintf(errors, "# Colonne: 1_errore posizione;   2_errore velocità\n#\n"); */
	
	K1 = malloc(4*sizeof(double));
	K2 = malloc(4*sizeof(double));
	
	J1 = malloc(4*sizeof(double));
	J2 = malloc(4*sizeof(double));

	t=0;	u=0;
	Err1=0;	Err2=0;
	fprintf(output, "%.4f\t%e\t%e\n", t, X1, X2);
	for(step=0; step<STEPS; step++)
	{
		temp1[step%2] = X1;
		temp2[step%2] = X2;
		ev_K(X1, X2, t, H, velocity, acceleration, K1, K2);
		X1 += deltaX(K1);
		X2 += deltaX(K2);
		t+=H;
		
		/* Stima dell'errore in funzione del passo reticolare H:
		 * confronto dell'evoluzione con 1 passo di ampiezza 2H con quella con
		 * 2 passi di ampiezza H
		 */
		if(((step%2)!=0)&&(STEPS%2==0))
		{
			ev_K(Z1, Z2, t, 2*H, velocity, acceleration, J1, J2);
			Z1 = temp1[0] + deltaX(J1);
			Z2 = temp2[0] + deltaX(J2);
			u+=2*H;
			Err1 += fabs(Z1 - X1);
			Err2 += fabs(Z2 - X2);
		}
		fprintf(output, "%.4f\t%e\t%e\n", t, X1, X2);
	}
	fprintf(errors, "%e\t%e\t%e\n", H, Err1, Err2);
	fclose(output);
	fclose(errors);
	
	exit(EXIT_SUCCESS);
}

