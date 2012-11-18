
/*******************************************************************************
 *
 *		File "integrali.c"
 *
 * Programma per il calcolo di integrali definiti con metodi deterministici.
 * I metodi utilizzati (definiti nel file "../modules/integration/methods.c")
 * sono in particolare:
 * _ metodo dei trapezi;
 * _ metodo delle parabole di Simpson;
 * _ metodo delle quadrature gaussiane (con polinomio di Legendre di grado 5).
 *
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integration.h"
#include "extras.h"

/*
 * Funzioni integrande
 */
double func_1 (double x)
{
	return log(1+x);
}

double func_2 (double x)
{
	return pow(x,9) - pow(x,7) + 3;
}


int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		printf("\nNumero insufficiente di parametri! Inserire:\n");
		printf("\t1_numero minimo di sottointervalli;\n");
		printf("\t2_numero massimo di sottointervalli;\n");
		printf("\t3_fattore di scala logaritmica.\n\n");
		exit(EXIT_FAILURE);
	}

	int n;
	int N_min = atoi(argv[1]);
	int N_max = atoi(argv[2]);
	int scale = atoi(argv[3]);

	/* Variabili per stime ed errori degli integrali */
	double I_trap, I_simp, I_gauss;
	double E_trap, E_simp, E_gauss;

	/* Risultati noti degli integrali (calcolati analiticamente) */
	double I1 = -1+log(27.0/4.0);
	double I2 = 73.425;

	/* Massimi dei moduli delle derivate di ordine 2, 4 e 10 delle funzioni
	 * integrande nell'intervallo di integrazione [1,2]
	 * (per la stima a priori dell'errore)
	 */
	double f1_2, f1_4, f1_10;
	double f2_2, f2_4, f2_10;	

	/* Puntatori alle funzioni integrande */
	double (*f1)(double);
	double (*f2)(double);
	f1 = func_1;
	f2 = func_2;

	f1_2  = 0.25;
	f1_4  = 0.375;
	f1_10 = 354.375;
	f2_2  = 7872.0;
	f2_4  = 90048.0;
	f2_10 = 0.0;

	/*
	 * File di output
	 */
	FILE *log_trap, *log_simp, *log_gauss;
		log_trap = fopen("IM/log_trap.dat","w");
		log_simp = fopen("IM/log_simp.dat","w");
		log_gauss = fopen("IM/log_gauss.dat","w");
	FILE *poly_trap, *poly_simp, *poly_gauss;
		poly_trap = fopen("IM/poly_trap.dat","w");
		poly_simp = fopen("IM/poly_simp.dat","w");
		poly_gauss = fopen("IM/poly_gauss.dat","w");
	
	/*
	 * Calcolo e stampa su file dei risultati e degli errori
	 * (a priori e reali)
	 */
	for(n = N_min; n < N_max + 1; n *= scale)
	{
		/*
		 * risultati per il logaritmo
		 */
		I_trap	= integration(f1, trapezoid, 1, 2, n);
		I_simp	= integration(f1, simpson, 1, 2, n);
		I_gauss	= integration(f1, gaussianquad_5, 1, 2, n);
		// formule per gli errori
		E_trap	= pow(1/(double)n,2)/12.0*f1_2;
		E_simp	= pow(1/(double)n,4)/90.0*f1_4;
		E_gauss	= pow(1/(double)n,10)/(double)factorial(10)*f1_10;
		
		fprintf(log_trap, "%d\t%.9f\t%e\t%e\n", n, I_trap, fabs(I_trap - I1), E_trap);
		fprintf(log_simp, "%d\t%.9f\t%e\t%e\n", n, I_simp, fabs(I_simp - I1), E_simp);
		fprintf(log_gauss, "%d\t%.9f\t%e\t%e\n", n, I_gauss, fabs(I_gauss - I1), E_gauss);
		
		/*
		 * risultati per il polinomio
		 */
		I_trap	= integration(f2, trapezoid, 1, 2, n);
		I_simp	= integration(f2, simpson, 1, 2, n);
		I_gauss	= integration(f2, gaussianquad_5, 1, 2, n);
		// formule per gli errori
		E_trap	= pow(1/(double)n,2)/12.0*f2_2;
		E_simp	= pow(1/(double)n,4)/90.0*f2_4;
		E_gauss	= pow(1/(double)n,10)/(double)factorial(10)*f2_10;

		fprintf(poly_trap, "%d\t%.9f\t%e\t%e\n", n, I_trap, fabs(I_trap - I2), E_trap);
		fprintf(poly_simp, "%d\t%.9f\t%e\t%e\n", n, I_simp, fabs(I_simp - I2), E_simp);
		fprintf(poly_gauss, "%d\t%.9f\t%e\t%e\n", n, I_gauss, fabs(I_gauss - I2), E_gauss);
	}
	fclose(log_trap);
	fclose(log_simp);
	fclose(log_gauss);
	fclose(poly_trap);
	fclose(poly_simp);
	fclose(poly_gauss);
	
	exit(EXIT_SUCCESS);
}
