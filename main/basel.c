/*******************************************************************************
 *
 *		File "basel.c"
 *
 * Programma per il calcolo numerico della serie di 1/n^2 in approssimazione di
 * somme parziali N-esime.
 * Il calcolo viene effettuato (con numeri sia in singola che in doppia
 * precisione) in due ordini di sommatoria opposti:
 * _ n crescente;
 * _ n decrescente.
 *
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 * Somma con indici crescenti con numeri double
 */
double dSUM_C (int N_max)
{	
	double S = 0;
	int i;
	for(i = 1; i < N_max+1; i++)
		S += 1/((double)i*(double)i);
		
	return S;
} 


/*
 * Somma con indici decrescenti con numeri double
 */
double dSUM_D (int N_max)
{		
	double S = 0;
	int i;
	for(i = N_max; i > 0; i--)
		S += 1/((double)i*(double)i);
	
	return S;
}


/*
 * Somma con indici crescenti con numeri float
 */
float fSUM_C (int N_max)
{	
	float S = 0;
	int i;
	for(i = 1; i < N_max+1; i++)
		S += 1/((float)i*(double)i);
	
	return S;
} 


/*
 * Somma con indici decrescenti con numeri float
 */
float fSUM_D (int N_max){
		
	float S = 0;
	int i;
	for(i = N_max; i > 0; i--)
		S += 1/((float)i*(double)i);
	
	return S;
}


int main(int argc, char * argv[])
{
	if(argc < 4)
	{
		printf("\nNumero insufficiente di parametri! Inserire:\n");
		printf("\t1_minimo indice di somma parziale;\n");
		printf("\t2_massimo indice di somma parziale;\n");
		printf("\t3_fattore di scala logaritmica.\n\n");
		exit(EXIT_FAILURE);
	}
	
	/*
	 * Variabili utili:
	 * _ N -> indice di somma parziale;
	 * _ N_min, N_max -> rispettivamente minimo e massimo indice di somma
	 *		parziale da calcolare;
	 * _ scale -> fattore di scala logaritmica con cui scegliere N;
	 * _ dTemp, fTemp -> arrays ausiliari.
	 */
	int N, N_min, N_max, scale;
	double dTemp[2];
	float fTemp[2];

	double PI;
	PI = 4.0*atan(1.0);

	/*
	 * File di output
	 */
	FILE *f_double, *f_float;
	f_double = fopen("serie/dati_double.dat", "w");
	f_float = fopen("serie/dati_float.dat","w");
	
	N_min = atoi(argv[1]);
	N_max = atoi(argv[2]);
	scale = atoi(argv[3]);
	
	/*
	 * Calcolo e stampa su file dei risultati e degli errori
	 * con i diversi ordini di sommatoria e al variare dell'indice N
	 */
	for(N = N_min; N < N_max; N *= scale)
	{
		dTemp[0] = dSUM_C(N);
		dTemp[1] = dSUM_D(N);
		fprintf(f_double, "%d\t%.9f\t%e\t%.9f\t%e\n", N, dTemp[0], PI*PI/6 - dTemp[0], dTemp[1], PI*PI/6 - dTemp[1]);
		
		fTemp[0] = fSUM_C(N);
		fTemp[1] = fSUM_D(N);
		fprintf(f_float, "%d\t%.9f\t%e\t%.9f\t%e\n", N, fTemp[0], PI*PI/6 - fTemp[0], fTemp[1], PI*PI/6 - fTemp[1]);
	}
	
	fclose(f_double);
	fclose(f_float);
	
	return(EXIT_SUCCESS);
}
