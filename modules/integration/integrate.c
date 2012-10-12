/*******************************************************************************
 * 
 *		File "integrate.c"
 *
 * Routine per calcolo di integrali definiti con metodi deterministici.
 * I metodi utilizzabili sono definiti nel file "methods.c".
 *
 ******************************************************************************/


#define INTEGRATE_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integration.h"

#define BINW 64

/* Routine di integrazione. Gli argomenti sono:
 * _ double (*f)(double) -> puntatore alla funzione integranda;
 * _ double (*method)(double (*g)(double), double, double) -> puntatore alla
 *		funzione per il metodo di integrazione;
 * _ double a, double b -> estremi di integrazione;
 * _ int N -> numero di sottointervalli per integrazione composita.
 */
double integration (double (*f)(double), double (*method)(double (*g)(double), double, double), double a, double b, int N)
{
	if(a == b)
		return 0;

	int i;
	double S = 0;
	double step;
	step = (b - a)/N;
	
	/* Per diminuire errori di somma (per N troppo grande) si binnano -se
	 * possibile- i sottointervalli.
	 */
	if((N%BINW)==0)
	{
		int m = N/BINW;
		for(i=0; i<BINW; i++)
			S += integration(f, method, a+i*m*step, a+(i+1)*m*step, m);
	}
	else
		for(i = 0; i< N; i++)
			S += method(f, a + i*step, a + step*(i + 1));
		
	return S;

}

