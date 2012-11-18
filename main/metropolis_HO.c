
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "random.h"
#include "simulations.h"
#include "cluster.h"

#include "extras.h"

#define delta	0.7		/* hypercube side */
#define NTH		200		/* number of sweeps needed for thermalization */
#define SWEEPS	1000000	/* number of "effective" sweeps */
#define NSKIP	100		/* number of sweeps skipped in order to avoid correlation */
#define TMAX	50		/* maximum time delay for autocorrelation */


/* Parameters for the trial wave functions */
double sigma;

/* Trial wave function and its transformed by hamiltonian operator */
double trial_WF_gauss (double x);
double H_trial_WF_gauss (double x);
double trial_WF_agnesi (double x);
double H_trial_WF_agnesi (double x);

/* pointers to functions corresponding to wave function and its transform respectively */
double (*WF)(double x);
double (*H_WF)(double x);

/* (Non normalized) probability distribution */
double probability (double *x);



int main (int argc, char *argv[])
{
	int sweep, choice, i;
	double x[1];
	/* double Sum;  Summation for the binning method */
	
	int N = SWEEPS/NSKIP;
	cluster exp_value;
		cluster_init(&exp_value, N);
	double *autocorr;
		autocorr = malloc(SWEEPS*sizeof(double));
	
	char *therm_name, *integral_name, *autocorr_name;
		therm_name		= malloc(100*sizeof(char));
		integral_name	= malloc(100*sizeof(char));
		autocorr_name	= malloc(100*sizeof(char));
	
	printf("\nChoose type of trial function:\n 1 _ gaussian,\t2 _ Agnesi's witch\t\t");
	scanf("%d", &choice);
	switch(choice)
	{		
		case 1:
			sprintf(therm_name, "MA_HO/gaussian_thermalization.dat");
			sprintf(integral_name, "MA_HO/gaussian_expectationvalues.dat");
			sprintf(autocorr_name, "MA_HO/gaussian_autocorrelation.dat");
			WF		= trial_WF_gauss;
			H_WF	= H_trial_WF_gauss;
			break;
		
		case 2:
			sprintf(therm_name, "MA_HO/witch_thermalization.dat");
			sprintf(integral_name, "MA_HO/witch_expectationvalues.dat");
			sprintf(autocorr_name, "MA_HO/witch_autocorrelation.dat");
			WF		= trial_WF_agnesi;
			H_WF	= H_trial_WF_agnesi;
			break;

		default:
			sprintf(therm_name, "MA_HO/gaussian_thermalization.dat");
			sprintf(integral_name, "MA_HO/gaussian_expectationvalues.dat");
			sprintf(autocorr_name, "MA_HO/gaussian_autocorrelation.dat");
			WF		= trial_WF_gauss;
			H_WF	= H_trial_WF_gauss;
	}
	
	FILE *therm_file, *integral_file, *autocorr_file;
		//therm_file		= fopen(therm_name, "w");
		integral_file	= fopen(integral_name, "w");
		autocorr_file	= fopen(autocorr_name, "w");

	printf("\nData will be saved in the following files:\n");
	printf("\t\"%s\"\t-> estimations of energy;\n", integral_name);
	printf("\t\"%s\"\t-> values of x at each sweep of the metropolis;\n", therm_name);
	printf("\t\"%s\"\t-> autocorrelation of data.\n", autocorr_name);
	
	/* For each value of the parameter of the trial wave function
	 * expectation value of the local energy is estimated
	 **/
	printf("Parameter of the trial wave function\n");
	for(sigma=0.5; sigma<1.6; sigma+=0.033)
	{	
		printf("\t%lf\n", sigma);
		
		/* Initialization of the randomizer */
		srand(time(NULL));
		rlxd_init(1,rand());
		
		/* first point is picked randomly ("hot" initialization) */
		hot_init(x,1);

		/* Process is left free for a certain number NTH of sweeps:
		 * no data are used for estimating the integral
		 **/
		for(sweep=0; sweep<NTH; sweep++)
		{
			metropolis(probability, x, 1, delta);
			/*if((sigma-0.5)/0.033 == 0)
				fprintf(therm_file, "%d\t%.10e\n", sweep, x[0]);*/
		}
		
		/* From now on data are collected and used for the estimation */
		i=0;
		for(sweep=0; sweep<SWEEPS; sweep++)
		{
			metropolis(probability, x, 1, delta);
			
			/* All data corresponding to each sweep of the metropolis algorithm
			* (apart from those needed for thermalization) are stored in the
			* autocorr array (in order to compute the autocorrelation).
			**/
			
			/*if((sigma-0.5)/0.033 == 0)
			{
				autocorr[sweep] = x[0];
				fprintf(therm_file, "%d\t%.10e\n", sweep+NTH, x[0]);
			}*/
			
			/* In order to get independent samples we take into account just
			 * the configuration every NSKIP sweeps; data are stored in an array
			 * contained in a jackknife structure.
			 **/
			if((sweep+1)%NSKIP==0)
			{
				exp_value.Vec[i] = H_WF(x[0])/WF(x[0]);
				i++;
			}
		}
				
		/* Expectation value and variance of the integral are computed
		 * with the jackknife re-sampling method
		 **/
		clusterJK(&exp_value);
		fprintf(integral_file, "%lf\t%.10e\t%.10e\n", sigma, exp_value.Mean, exp_value.Sigma);
	}
	
	for(i=0; i<TMAX+1; i++)
			fprintf(autocorr_file, "%d\t%e\n", i, autocorrelation(autocorr, i, N));
	
	//fclose(therm_file);
	fclose(autocorr_file);
	fclose(integral_file);
	free(autocorr);
	
	exit(EXIT_SUCCESS);
}



double trial_WF_gauss (double x)
{	return exp(-x*x/(2*sigma*sigma)); }


double H_trial_WF_gauss (double x)
{	return exp(-x*x/(2*sigma*sigma))*(sigma*sigma + x*x*(pow(sigma,4)-1))/2/pow(sigma,4);  }


double trial_WF_agnesi (double x)
{	return 1/(x*x + sigma*sigma);  }


double H_trial_WF_agnesi (double x)
{	return (2*sigma*sigma + x*x*((x*x+sigma*sigma)*(x*x+sigma*sigma) - 6))/(2*((x*x+sigma*sigma)*(x*x+sigma*sigma)*(x*x + sigma*sigma))); }


double probability (double *x)
{	return WF(x[0])*WF(x[0]);  }
