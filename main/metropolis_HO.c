
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "random.h"
#include "metropolis.h"
#include "cluster.h"

#define delta	2.0		/* hypercube side */
#define NTH		200		/* number of sweeps needed for thermalization */
#define SWEEPS	1000000	/* number of "effective" sweeps */
#define NSKIP	100		/* number of sweeps skipped in order to avoid correlation */
#define TMAX	20		/* maximum time delay for autocorrelation */


/* Parameters for the trial wave functions */
double sigma;
//double alpha;

/* Trial wave function and its transformed by hamiltonian operator */
double trial_WF (double x);
double H_trial_WF (double x);

/* (Non normalized) probability distribution */
double probability (double x);



int main (int argc, char *argv[])
{
	int sweep, choice, i;
	double x[1];
	double Sum, Sum1;
	
	int N = SWEEPS/NSKIP;
	cluster exp_value;
		cluster_init(&exp_value, N);
	
	char *therm_name, *integral_name, *autocorr_name;
		therm_name		= malloc(100*sizeof(char));
		integral_name	= malloc(100*sizeof(char));
		autocorr_name	= malloc(100*sizeof(char));
	
	printf("\nChoose type of trial function:\n 1 _ gaussian,\t2 _ Agnesi's witch\t\t");
	switch(choice)
	{
		case 1:
			sprintf(therm_name, "MA_HO/gaussian_thermalization.dat");
			sprintf(integral_name, "MA_HO/gaussian_expectationvalues.dat");
			sprintf(integral_name, "MA_HO/gaussian_autocorrelation.dat");
		
		case 2:
			sprintf(therm_name, "MA_HO/witch_thermalization.dat");
			sprintf(integral_name, "MA_HO/witch_expectationvalues.dat");
			sprintf(integral_name, "MA_HO/witch_autocorrelation.dat");
			
		default:
			sprintf(therm_name, "MA_HO/gaussian_thermalization.dat");
			sprintf(integral_name, "MA_HO/gaussian_expectationvalues.dat");
			sprintf(integral_name, "MA_HO/gaussian_autocorrelation.dat");
	}
	
	FILE *therm_file, *integral_file, *autocorr_file;
		therm_file		= fopen(therm_name, "w");
		integral_file	= fopen(integral_name, "w");
		autocorr_file	= fopen(autocorr_name, "w");
	
	for(sigma=0.5; sigma<1.6; sigma+=0.05)
	{
		i=0;
		srand(time(NULL));
		rlxd_init(1,rand());
		
		hot_init(x,1);
		
		for(sweep=0; sweep<NTH; sweep++)
		{
			metropolis(probability, x, 1, delta);
			if((sigma-0.5)/0.05 == 1)
				fprintf(therm_file, "%d\t%.10e\n", sweep, x[0]);
		}
		
		for(sweep=0; sweep<SWEEPS; sweep++)
		{
			metropolis(probability, x, 1, delta);
			if((sigma-0.5)/0.05 == 1)
				fprintf(therm_file, "%d\t%.10e\n", sweep+NTH, x[0]);
			
			/* HERE COME THE MAIN STEPS !!!!! */
			if((sweep+1)%NSKIP==0)
			{	
				exp_value.Vec[i] = H_trial_WF(x[0])/trial_WF(x[0]);
				i++;
			}
		}
		clusterJK(&exp_value);
		fprintf(integral_file, "%lf\t%.10e\t.10e\n", sigma, exp_value.Mean, exp_value.Sigma);
	}
	exit(EXIT_SUCCESS);
}



double trial_WF (double x)
{
	return exp(x*x/(2*sigma*sigma));
	//return 1/(x*x + alpha*alpha)
}


double H_trial_WF (double x)
{
	return exp(x*x/(2*sigma*sigma))*(sigma*sigma + x*x*(pow(sigma,4)-1))/2/pow(sigma,4);
	
	/* return (2*alpha*alpha + x*x*((x*x+alpha*alpha)*(x*x+alpha*alpha) - 6))/(2*((x*x+alpha*alpha)*(x*x+alpha*alpha)*(x*x+alpha*alpha))); */
}


double probability (double x)
{
	return trial_WF(x)*trial_WF(x);
}
