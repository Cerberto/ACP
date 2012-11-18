/*******************************************************************************
 * 
 * 		File "harmosc.c"
 * 
 * 
 * Routines necessarie per la simulazione numerica con algoritmo di Metropolis
 * dell'oscillatore armonico quantistico.
 * 
 * 
 * Le routines accessibili dall'esterno sono:
 * _ HOeAction -> azione euclidea dell'oscillatore armonico;
 * 
 * _ HOmetropolis -> algoritmo di scelta di Metropolis (restituisce la
 * 		variazione di azione conseguente la scelta);
 *
 * _ DeltaE -> calcolo del gap di energia (restituisce una struttura cluster
 * 		jackknife contenente media ed errore);
 * 
 * _ MatrixElementX -> calcolo dell'elemento di matrice dell'operatore posizione
 * 		tra autostati fondamentale e primo eccitato (restituisce una struttura
 * 		cluster jackknife contenente media ed errore).
 *
 * _ sqrt_jk -> calcolo della radice quadrata su un campione di dati campionati
 *		col metodo jackknife
 * 
 * 
 * Le routines private sono:
 * _ potential -> potenziale parabolico;
 * 
 * _ HOeLagrangian -> lagrangiana euclidea dell'oscillatore armonico;
 * 
 * _ deltaS -> variazione di azione tra due configurazioni con una sola
 * 		"componente" differente;
 * 
 ******************************************************************************/


#define HARMOSC_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "simulations.h"
#include "random.h"
#include "cluster.h"

#define M 1.0
#define OMEGA 1.0
#define DELTA 4.0


/* Potenziale parabolico */
static double potential(double x)
{
	return M*OMEGA*OMEGA*x*x/2;
}


/* Lagrangiana euclidea */
static double HOeLagrangian(double x1, double x2)
{
	return M/2*(x2 - x1)*(x2 - x1) + (potential(x1) + potential(x2))/2;
}


/* Azione euclidea */
double HOeAction(double * x, int dim)
{
	int j;
	double S=0;
	for(j=0; j<dim; j++)
		S += HOeLagrangian(x[(j+1)%dim], x[j]);
	
	return S;
}


/* Variazione di azione euclidea conseguente alla modifica di un elemento del
 * vettore x
 */
static double deltaS(double *x, int dim, double x_new, int i)
{
	return M*((x[i] - x_new)*(x[(i+1)%dim] + x[(i-1+dim)%dim]) +\
	 (x_new*x_new - x[i]*x[i])) + potential(x_new) - potential(x[i]);
}


/* Algoritmo di Metropolis */
double HOmetropolis (double* state, int state_dim)
{
	int i;
	double u[2];
	double x_new, temp1, temp2;
	double DS = 0;
	
	for(i=0; i<state_dim; i++)
	{
		ranlxd(u,2);
		x_new = state[i] + DELTA*(2*u[1] - 1);
		temp1 = deltaS(state,state_dim,x_new,i);
		temp2 = exp(-temp1);

		/* scelta nuova configurazione */
		if(temp2 >= u[0])
		{
			state[i] = x_new;
			DS += temp1;
		}
	}
	return DS;
}


void HOautocorrelation(double *x, int t, double *v, int dim, int steps)
{
	int i;
	double *temp1, *temp2, *temp3;
	temp1 = malloc(dim*sizeof(double));
	temp2 = malloc(dim*sizeof(double));
	temp3 = malloc(dim*sizeof(double));
	cold_init(temp1, dim);
	cold_init(temp2, dim);
	cold_init(temp3, dim);
	
	for(i=0; (i+dim*t)<(steps*dim); i++)
	{
		temp3[i%dim] += x[i]*x[i+t*dim]/(double)(steps - t);
		temp1[i%dim] += x[i]/(double)(steps - t);
		temp2[i%dim] += x[i]*x[i]/(double)(steps - t);
	}
	for(i=0; i<dim; i++)
		v[i] = (temp3[i] - temp1[i]*temp1[i])/(temp2[i] - temp1[i]*temp1[i]);
	
	free(temp1);
	free(temp2);
	free(temp3);
}



/* Gap di energia */
cluster DeltaE(cluster *A, cluster *B, cluster *C)
{
	if(((A->Dim)!=(B->Dim))||((A->Dim)!=(C->Dim)))
	{
		printf("\nClusters ad argomento della funzione DeltaE ");
		printf("contengono array di diversa lunghezza!\n\n");
		exit(EXIT_FAILURE);
	}
	
	int i;
	double temp = 0;
	int dim = (A->Dim);
	cluster result;
	cluster_init(&result, dim);
	result.Mean	= acosh((A->Mean + C->Mean)/2.0/(B->Mean));
	for(i=0; i<dim; i++)
	{
		result.Vec[i] = acosh((A->Vec[i] + C->Vec[i])/2.0/(B->Vec[i]));
		temp += ((double)(dim - 1)/(double)dim)*(result.Vec[i] - result.Mean)*(result.Vec[i] - result.Mean);
	}
	result.Sigma = temp;
	
	return result;
}


/* Elemento di matrice operatore posizione */
cluster MatrixElementX(cluster *DE, cluster *Corr, int t, int N)
{
	if((DE->Dim)!=(Corr->Dim))
	{
		printf("\nClusters ad argomento della funzione MatrixElementX ");
		printf("contengono array di diversa lunghezza!\n\n");
		exit(EXIT_FAILURE);
	}
	
	int i;
	double temp = 0;
	int dim = DE->Dim;
	cluster result;
	cluster_init(&result,dim);
	result.Mean = (Corr->Mean)*exp(0.5*N*(DE->Mean))/cosh((0.5*N - t)*(DE->Mean));
	result.Mean = sqrt(result.Mean);
	for(i=0; i<dim; i++)
	{
		result.Vec[i] = (Corr->Vec[i])*exp(0.5*N*(DE->Vec[i]))/cosh((0.5*N - t)*(DE->Vec[i]));
		result.Vec[i] = sqrt(result.Vec[i]);
		temp += ((double)(dim - 1)/(double)dim)*(result.Vec[i] - result.Mean)*(result.Vec[i] - result.Mean);
	}
	result.Sigma = temp;
	
	return result;
}


/* Funzione radice quadrata: media e varianza su un campione */
cluster sqrt_jk (cluster *X)
{
/*
	int i;
	double temp = 0;
	int dim = X->Dim;
	cluster result;
	cluster_init(&result,dim);
	result.Mean = sqrt(X->Mean);
	for(i=0; i<dim; i++)
	{
		result.Vec[i] = sqrt(X->Vec[i]);
		temp += ((double)(dim - 1)/(double)dim)*(result.Vec[i] - result.Mean)*(result.Vec[i] - result.Mean);
	}
	result.Sigma = temp;
	
	return result;
*/
	return functionJK(sqrt, X);
}
