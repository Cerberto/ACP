
/*******************************************************************************
 *
 *		File "expo.c"
 *
 *	Generazione di numeri con distribuzione esponenziale.
 * Le funzioni accessibili dall'esterno sono:
 * _ symrootexpo, symrootexpo_dble -> generano numeri random (rispettivamente a
 *		singola e doppia precisione) con distribuzione proporzionale a
 *		sqrt(|x|)*e^(-|x|) e li assegnano agli elementi del vettore r passato ad
 *		argomento;
 * _ symrootexpodistr -> funzione esponenziale decrescente (corrispondente alla
 *		distribuzione dei numeri generati da symexpo e symexpo_double;
 *
 ******************************************************************************/

#define SYMROOTEXPO_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"

#define PI 3.141592653589793

void symrootexpo (float vf[], int n)
{
	int k;
	float temp_exp[1];
	float temp_gauss[1];
	double x1, x2, y;
	double choice[1];
	
	for(k=0; k<n; k++)
	{	
		ranlxd(choice, 1);
		gauss(temp_gauss,1);
		expo(temp_exp,1);
		x1 = (double)temp_gauss[0];
		x2 = (double)temp_exp[0];
		y = pow(x1,2) + x2;
		
		if (choice[0] > 0.5)
			vf[k] = (float)y;
		
		else
			vf[k] = -(float)y;
	}
}


void symrootexpo_dble (double vd[], int n)
{
	int k;
	double temp_exp[1];
	double temp_gauss[1];
	double x1, x2, y;
	
	for(k=0; k<n; k++)
	{
		gauss_dble(temp_gauss,1);
		expo_dble(temp_exp,1);
		x1 = temp_gauss[0];
		x2 = temp_exp[0];
		y = pow(x1,2) + x2;

		vd[k] = y;
	}
}


double symrootexpodistr (double x)
{
	return sqrt(x/PI)*exp(-x);
}

