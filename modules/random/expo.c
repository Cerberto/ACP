
/*******************************************************************************
 *
 *		File "expo.c"
 *
 *	Generazione di numeri con distribuzione esponenziale.
 * Le funzioni accessibili dall'esterno sono:
 * _ expo, expo_dble -> generano numeri random (rispettivamente a singola e
 *		doppia precisione) positivi con distribuzione esponenziale e li
 *		assegnano agli elementi del vettore r passato ad argomento;
 * _ expodistr -> funzione esponenziale decrescente (corrispondente alla
 *		distribuzione dei numeri generati da expo e expo_double;
 * _ wexpo, wexpo_dble -> generano numeri random (rispettivamente a singola
 *		e doppia precisione) positivi con distribuzione proporzionale a
 *		e^(-x/alpha) (dove alpha è passato ad argomento) e li assegnano agli
 *		elementi del vettore r passato ad argomento;
 * _ wexpodistr -> funzione corrispondente alla distribuzione dei numeri
 *		generati da wexpo e wexpo_dble.
 ******************************************************************************/


#define EXPO_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
#include "extras.h"

#define PI 3.141592653589793


/*
 * Distribuzione esponenziale
 */

void expo (float r[], int n)
{
	int k;
	float u[1];
	double x, y;
	

	for (k=0; k<n; k++)
	{
		ranlxs(u,1);
		x=(double)u[0];

		y = -log(1.0-x);
      
		r[k]=(float)y;
	}
}


void expo_dble (double r[], int n)
{
	int k;
	double u[1];
	double x, y;

	for (k=0; k<n; k++)
	{
		ranlxd(u,1);
		x = u[0];
      
		y = -log(1.0-x);
      
		r[k] = y;
	}
}


double expodistr (double x)
{
	return exp(-x);
}




/*
 * Distribuzioni esponenziali con larghezza definita dall'esterno
 */

void wexpo (float r[], int n, float alpha)
{
	int k;
	float u[1];
	double x, y;
	

	for (k=0; k<n; k++)
	{
		ranlxs(u,1);
		x=(double)u[0];

		y = -(double)alpha*log(1.0-x);
      
		r[k]=(float)y;
	}
}


void wexpo_dble (double r[], int n, double alpha)
{
	int k;
	double u[1];
	double x, y;

	for (k=0; k<n; k++)
	{
		ranlxd(u,1);
		x = u[0];
      
		y = -alpha*log(1.0-x);
      
		r[k] = y;
	}
}


double wexpodistr (double x, double alpha)
{
	if (alpha <= 0)
	{
		printf("\nNon-positive width parameter!!\n\n");
		exit(EXIT_FAILURE);
	}
	
	else
		return exp(-x/alpha)/alpha;
}
