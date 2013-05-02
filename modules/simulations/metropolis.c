

#define METROPOLIS_C

#include <stdio.h>
#include <stdlib.h>
#include "simulations.h"
#include "random.h"
#include "cluster.h"
#include "extras.h"

/* Cold initialization */
void cold_init(double *v, int dim)
{
	int i;
	for(i=0; i<dim; i++)
		v[i] = 0;
}


/* Hot initialization */
void hot_init(double *v, int dim)
{
	int i;
	double *temp;
	temp = malloc(dim*sizeof(double));
	ranlxd(temp,dim);
	for(i=0; i<dim; i++)
		v[i] = 10.0*(2*temp[i] - 1);
	
	free(temp);
}



/* Autocorrelation of data in an array */
double autocorrelation(double *x, int t, int dim)
{
	int i;
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;

	for(i=0; i<(dim-t); i++)
	{
		temp3 += x[i]*x[i+t]/(double)(dim - t);
		temp1 += x[i]/(double)(dim - t);
		temp2 += x[i]*x[i]/(double)(dim - t);
	}

	return (temp3 - temp1*temp1)/(temp2 - temp1*temp1);
}


/* Routine that executes a sweep of the Metropolis algorithm */
void metropolis (double (*P)(double *), double *state, int state_dim, double delta)
{
	int i;
	double swap, x_new, acceptance;
	double u[2];

	for(i=0; i<state_dim; i++)
	{
		ranlxd(u,2);
		x_new = state[i] + delta*(u[0] - 0.5);
		swap = state[i];
		state[i] = x_new;
			acceptance = P(state);
		state[i] = swap;
			acceptance /= P(state);
		
		if(acceptance >= u[1])
			state[i] = x_new;
	}
}
