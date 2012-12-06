
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "random.h"
#include "simulations.h"
#include "cluster.h"
#include "extras.h"

#define NTH		200		/* number of sweeps needed for thermalization */
#define SWEEPS	100000	/* number of "effective" sweeps */
#define NSKIP	100		/* number of sweeps skipped in order to avoid correlation */
#define TMAX	150		/* maximum time delay for autocorrelation */

/* hypercube side */
double delta = 0.2;

/* Parameters for the trial wave functions */
double sigma;

/* Trial wave function and its transformed by hamiltonian operator */
double trial_WF (double *x, int dim);
double H_trial_WF (double *x, int dim);

/* pointers to functions corresponding to wave function and its transform respectively */
double (*WF)(double *x, int dim);
double (*H_WF)(double *x, int dim);

/* (Non normalized) probability distribution */
double probability (double *x, int dim);



int main (int argc, char *argv[])
{	
	int sweep, choice, i;
	double x[1];
	double sigma_min, sigma_max;
	/* double Sum;  Summation for the binning method */
	
	int N = SWEEPS/NSKIP;
	cluster exp_value;
		cluster_init(&exp_value, N);
	double *autocorr;
		autocorr = malloc(SWEEPS*sizeof(double));
	
	char *type, *distrib_name, *therm_name, *integral_name, *autocorr_name;
		type			= malloc(10*sizeof(char));
		distrib_name	= malloc(100*sizeof(char));
		therm_name		= malloc(100*sizeof(char));
		integral_name	= malloc(100*sizeof(char));
		autocorr_name	= malloc(100*sizeof(char));
	
	FILE *distrib_file, *therm_file, *integral_file, *autocorr_file;
		distrib_file	= fopen(distrib_name, "w");
		therm_file		= fopen(therm_name, "w");
		integral_file	= fopen(integral_name, "w");
		autocorr_file	= fopen(autocorr_name, "w");

	printf("\nData will be saved in the following files:\n");
	printf("\t\"%s\"\t-> estimations of energy;\n", integral_name);
	printf("\t\"%s\"\t-> values of x at each sweep of the metropolis;\n", therm_name);
	printf("\t\"%s\"\t-> values of x when system is assumed to be \"thermalized\";\n", distrib_name);
	printf("\t\"%s\"\t-> autocorrelation of data.\n\n", autocorr_name);
	
	/* For each value of the parameter of the trial wave function
	 * expectation value of the local energy is estimated
	 **/
	printf("Parameter of the trial wave function\n");
	for(sigma=sigma_min; sigma<sigma_max; sigma+=0.033){	
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
			if((int)((sigma-sigma_min)/0.033) == 5)	
				fprintf(therm_file, "%d\t%.10e\n", sweep, x[0]);
		}
		if((int)((sigma-sigma_min)/0.033) == 5)
			printf("\tHere I print something on file, to check the distribution.\n");
		
		/* From now on data are collected and used for the estimation */
		i=0;
		for(sweep=0; sweep<SWEEPS; sweep++){
			metropolis(probability, x, 1, delta);
			
			/* All data corresponding to each sweep of the metropolis algorithm
			* (apart from those needed for thermalization) are stored in the
			* autocorr array (in order to compute the autocorrelation), just
			* for a particular value of the parameter sigma
			**/
			
			if((int)((sigma-sigma_min)/0.033) == 5){
				autocorr[sweep] = x[0];
				fprintf(distrib_file, "%.10e\n", x[0]);
				fprintf(therm_file, "%d\t%.10e\n", sweep+NTH, x[0]);
			}
			
			/* In order to get independent samples we take into account just
			 * the configuration every NSKIP sweeps; data are stored in an array
			 * contained in a jackknife structure.
			 **/
			if((sweep+1)%NSKIP==0){
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
	
	fclose(therm_file);
	fclose(distrib_file);
	fclose(autocorr_file);
	fclose(integral_file);
	free(autocorr);
	
	exit(EXIT_SUCCESS);
}



double trial_WF (double *x, int dim){
	return 1;
}


/* Just the laplacian of the wave function at the moment.... */
double H_trial_WF (double *x, int dim){
	
	int i,j,k;
	double *distances;
		distances = malloc((NPAR*(NPAR-1)/2)*sizeof(double));
	double temp;
	double S, sum1[3], sum2, sum3;
	
	/* Filling array containing mutual distances;
	 * distance between i-th and j-th particles (i and j run between 0 and (N-1))
	 * is in the position [j(j-1)/2 +i].
	 **/
	for(j=1; j<NPAR; j++){
		for(i=0; i<j; i++){
			S = 0;
			for(k=0; k<3; k++)
				S += (x[3*i+k] - x[3*j+k])*(x[3*i+k] - x[3*j+k]);
			distances[j*(j-1)/2 + i] = sqrt(S);
		}
	}
	
	/* Variables containing accumulating sums are initially set to 0 */
	S = 0;
	for(k=0; k<3; k++)
		sum1[k]=0;
	sum2 = 0;
	sum3 = 0;
	
	/* kinetic terms */
	for(j=0; i<NPAR; j++){
		for(i=0; i<NPAR; i++){
			if(i==l) continue;
			
			temp	= pow(distances[j(j-1)/2 +i], 7);
			sum1[k]	+= (x[3*i +k] - x[3*j +k])/temp;
			sum2	+= 1/temp;
		}
		for(k=0; k<3; k++)
			sum3 += sum1[k]*sum1[k];
		S += -1/2/M*(25.0*pow(B,10)/4.0*sum3 - 10*pow(B,5)*sum2)
	}
	
	/* interaction potential */
	for(j=1; j<NPAR; j++){
		for(i=0; i<j; i++){
			temp = pow(SIGMA/distances[j(j-1)/2 +i], 6);
			sum3 += (temp - temp*temp);
		}
	}
	
	/* external (harmonic) potential */
	for(i=0; i<NPAR; i++){
		temp = 0;
		for(k=0; k<3; k++)
			temp += x[3*i +k]*x[3*i +k];
		sum2 += temp;
	}
	
	S += 0.5*M*OMEGA*OMEGA*sum2 + 4*EPS*sum3;
	
	free(distances);
	
	return S;
}


double probability (double *x){
	return WF(x[0])*WF(x[0]);
}
