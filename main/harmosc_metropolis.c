
/*******************************************************************************
 * 
 * 		File "harmosc_metropolis.c"
 * 
 * Programma di simulazione numerica con algoritmo di Metropolis
 * dell'oscillatore armonico quantistico.
 *
 * Le funzioni esterne utilizzate sono definite in:
 * _ ../modules/random (generazione numeri casuali)
 *			ranlxd.c
 * _ ../modules/simulations
 *			harmosc.c   (simulazione oscillatore armonico)
 * _ ../modules/cluster  (routines per ricampionamento jackknife)
 *			jackknife.c
 * 
 ******************************************************************************/


#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "random.h"
#include "simulations.h"
#include "cluster.h"


/* Numero passi reticolo */
#define N 32

/* Massimo valore di Dt per cui calcolare l'autocorrelazione */
#define DMAX 4

/* Numero di correlatori utilizzati per calcolo di DeltaE */
#define NCL 5

/* Numero di sweeps Metropolis per termalizzazione */
#define NTH 200

/* Dimensione dei bin */
#define DBIN 100

/* Valore massimo di t, variabile dell'autocorrelazione */
#define TMAX 20

/* Numero di sweeps Metropolis "effettivi", i.e., utilizzati per il calcolo */
int STEPS = 100000;


int main (int argc, char * argv[])
{
	if(argc == 1)
	{
		printf("\nNumero di parametri insufficiente!");
		printf("\nImpostare tipo di inzializzazione:\n");
		printf("\t1 -> a freddo;\n");
		printf("\t2 -> a caldo\n\n");
		exit(EXIT_FAILURE);
	}

	if(argc > 2)
	{
		STEPS = atoi(argv[2]);
		if(STEPS%DBIN != 0)
		{
			printf("\n# sweeps non è multiplo intero della lunghezza dei bins!\n\n");
			exit(EXIT_FAILURE);
		}
	}

	int Nbins = STEPS/DBIN;
	int i, step, Dt, bin, t;
	int choice = 0;
	double action, Sum, Sum1, Err, Err1;
	
	/* Stringhe per il nome del file di output */
	char *file_autocorr, *file_action;
	file_autocorr	= malloc(100*sizeof(char));
	file_action		= malloc(100*sizeof(char));

	/*
	 * Vettori utili:
	 * _ state -> contiene la posizione dell'oscillatore in ogni passo del 
	 * 		reticolo (N passi totali);
	 * 
	 * _ corrDtstep -> contiene i valori stimati a ciascun passo del Metropolis 
	 * 		di <x_i*x_{i+Dt}>;
	 * 
	 * _ Vtemp, Vtemp1 -> vettori ausiliari per costruire i correlatori di x e
	 *		di x^2 con il metodo del binning;
	 **/
	double *state, *corrDtstep, *Vtemp, *Vtemp1, *autocorr;

	/*
	 * Clusters jackknife utili:
	 * _ clusterDt -> N cluster jk che contengono i correlatori di x e gli
	 *		errori;
	 * 
	 * _ clusterSQDt -> cluster jk che contengono i correlatori di x^2 e gli
	 *		errori;
	 *
	 * _ DEtemp e EMtemp -> clusters jk ausiliari per il calcolo del gap di
	 *		energia e dell'elemento di matrice <0|x|1>
	 **/
	cluster *clusterDt, *clusterSQDt;
	cluster DEtemp, EMtemp;

	/* Eliminazione dei vecchi dati e creazione cartelle per i nuovi dati
	 * dell'autocorrelazione
	 */
	system("rm -r harmosc/autocorrelation");
	system("mkdir harmosc/autocorrelation");

	/* 
	 * File di output utili:
	 * _ out_action -> azione dell'oscillatore armonico in funzione dello sweep
	 * 		del Metropolis;
	 * 
	 * _ out_corr -> correlatori <x_l x_k> ed errori;
	 * 
	 * _ out_deltaE -> valori di deltaE (mediato sugli Nbins bin) in funzione di
	 * 		Dt (variabile della correlazione);
	 * _ out_deltaE_errors -> errori sul calcolo di DeltaE;
	 * 
	 * _ out_EM -> valore elemento di matrice di <0|x|1>;
	 * _ out_EM_errors -> errori sul calcolo di <0|x|1>;
	 *
	 * _ out_x2gs -> risultati ed errori di <x^2> sullo stato
	 *		fondamentale;
	 * 
	 * _ out_corrsq -> correlatori <x^2_l x^2_k> ed errori;
	 * 
	 * _ out_autocorr -> files labellati dal valore di Dt in cui viene salvata
	 * 		la funzione di autocorrelazione per <x_l x_{l+Dt}>
	 */
	FILE *out_action, *out_corr, *out_deltaE, *out_EM, *out_deltaE_errors,\
		*out_EM_errors, *out_x2gs, *out_corrsq;
	FILE *out_autocorr[DMAX+1];
	out_corr	= fopen("harmosc/corr_results.dat","w");
	out_deltaE	= fopen("harmosc/deltaE_results.dat","a");
	out_EM		= fopen("harmosc/matrixelement_results.dat","a");
	out_deltaE_errors	= fopen("harmosc/deltaE_errors.dat","a");
	out_EM_errors		= fopen("harmosc/matrixelement_errors.dat","a");
	out_x2gs	= fopen("harmosc/Egs_results.dat", "a");
	out_corrsq	= fopen("harmosc/corrSQ_results.dat", "a");
	for(i=0; i<DMAX+1; i++)
	{
		sprintf(file_autocorr, "harmosc/autocorrelation/autocorr_%d", i);
		out_autocorr[i] = fopen(file_autocorr, "w");
	}

	/* Allocazione vettori e cluster jk utili */
	state 		= malloc(N*sizeof(double));
	corrDtstep	= malloc((N*STEPS)*sizeof(double));
	Vtemp		= malloc(N*sizeof(double));
	Vtemp1		= malloc(N*sizeof(double));
	autocorr	= malloc(N*sizeof(double));
	clusterDt	= malloc(N*sizeof(cluster));
	clusterSQDt	= malloc(N*sizeof(cluster));

	/* Inizializzazione strutture cluster jackknife */
	for(Dt=0; Dt<N; Dt++)
	{
		cluster_init(clusterDt+Dt,Nbins);
		cluster_init(clusterSQDt+Dt,Nbins);
	}
	cluster_init(&DEtemp,Nbins);
	cluster_init(&EMtemp,Nbins);

	srand(time(NULL));
	rlxd_init(1,rand());

	/* Scelta di inizializzazione "a freddo" o "a caldo" dello stato iniziale */
	choice = atoi(argv[1]);
	switch(choice)
	{
		case 1:
			cold_init(state,N);
			sprintf(file_action, "harmosc/action_coldinit.dat");
			break;

		case 2:
			hot_init(state,N);
			sprintf(file_action, "harmosc/action_hotinit.dat");
			break;
	}
	out_action	= fopen(file_action,"w");

	cold_init(Vtemp,N);

	fprintf(out_action,\
		"#\n# Azione euclidea ad ogni sweep dell'algoritmo Metropolis\n");
	fprintf(out_action,"# Le colonne sono:\n");
	fprintf(out_action,"# Sweep\t azione\n#\n");


	/* Valuto l'azione del sistema nella configurazione iniziale
	 * (dipenderà dal tipo di inizializzazione di state scelta)
	 */
	action = HOeAction(state,N);
	fprintf(out_action,"%d\t%e\n", -1, action);

	/* Finché non si è raggiunto il tempo di termalizzazione NTH faccio evolvere
	 * il sistema con il Metropolis
	 */
	for(step=0; step<NTH; step++)
	{
		action += HOmetropolis(state,N);
		fprintf(out_action,"%d\t%e\n", step, action);
	}

	bin = 0;
	/* Raggiunto il tempo di termalizzazione, si fa evolvere il sistema per un
	 * numero STEPS di sweeps del Metropolis utilizzando gli stati ottenuti
	 * per il calcolo della correlazione
	 */
	for(step=NTH; step<STEPS+NTH; step++)
	{
		action += HOmetropolis(state,N);
		fprintf(out_action,"%d\t%e\n", step, action);

		for(Dt=0; Dt<N; Dt++)
		{
			Sum = 0;
			Sum1 = 0;
			/* media sul vettore di reticolo (somme su i) di x_i*x_{i+K} e di
			 * x^2_i*x^2_{i+K}
			 */
			for(i=0; i<N; i++)
			{
				Sum  += state[i]*state[(i+Dt)%N];
				Sum1 += state[i]*state[(i+Dt)%N]*state[i]*state[(i+Dt)%N];
			}
			Vtemp[Dt]  += Sum/((double)N);
			Vtemp1[Dt] += Sum1/((double)N);
		}

		/* Salvataggio della media sul bin nell vettore del cluster jackknife */
		if((step-NTH+1)%DBIN == 0)
		{
			for(Dt=0; Dt<N; Dt++)
			{
				clusterDt[Dt].Vec[bin] = Vtemp[Dt]/(double)DBIN;
				Vtemp[Dt] = 0;
				clusterSQDt[Dt].Vec[bin] = Vtemp1[Dt]/(double)DBIN;
				Vtemp1[Dt] = 0;
			}
			bin++;
		}
	}

	Sum = 0; Err = 0;
	/* Calcolo e stampa di media e deviazione standard della media dei 
	 * correlatori.
	 */
	for(Dt=0; Dt<N; Dt++)
	{
		clusterJK(clusterDt+Dt);
		fprintf(out_corr,"%d\t%e\t%e\n", Dt, clusterDt[Dt].Mean, \
			sqrt(clusterDt[Dt].Sigma));
		clusterJK(clusterSQDt+Dt);		
		fprintf(out_corrsq, "%d\t%e\t%e\n", Dt, clusterSQDt[Dt].Mean, \
			sqrt(clusterSQDt[Dt].Sigma));
		
		/* Calcolo l'elemento di matrice di x^2 sul ground state */
		EMtemp = sqrt_jk(clusterSQDt+Dt);
		Sum		+= (EMtemp.Mean)/(EMtemp.Sigma);
		Err		+= 1.0/(EMtemp.Sigma);
	}
	fprintf(out_x2gs, "%e\t%e\n",Sum/Err, sqrt(Err/((double)Nbins)));
	fprintf(out_x2gs, "%e\n", Sum/Err);

	/* Calcolo del gap di energia e degli elementi di matrice con varianze.
	 * Si prendono in considerazione soltanto i correlatori per piccoli
	 * valodi di Dt: i valori centrali non hanno un andamento regolare
	 * (si veda dal plot dei correlatori)
	 */
	for(Dt=2; Dt<NCL; Dt++)
	{
		Sum = 0; Sum1 = 0;
		Err = 0; Err1 = 0;
		
		DEtemp = DeltaE(clusterDt+(Dt-1), clusterDt+Dt, clusterDt+(Dt+1));
		EMtemp = MatrixElementX(&DEtemp, clusterDt+Dt, Dt, N);
			
		/* Si effettua la media pesata di DeltaE per i valori di Dt
		 * presi in considerazione
		 */
		Sum		+= (DEtemp.Mean)/(DEtemp.Sigma);
		Sum1	+= (EMtemp.Mean)/(EMtemp.Sigma);

		Err		+= 1.0/(DEtemp.Sigma);
		Err1	+= 1.0/(EMtemp.Sigma);
	}

	fprintf(out_deltaE_errors,	"%d\t%e\n", STEPS, sqrt(Err/((double)Nbins)));
	fprintf(out_EM_errors,		"%d\t%e\n", STEPS, sqrt(Err1/((double)Nbins)));

	fprintf(out_deltaE,	"%e\n",Sum/Err);
	fprintf(out_EM,		"%e\n",Sum1/Err1);

	/*
	 * Calcolo e stampa delle autocorrelazioni per i correlatori delle x per i
	 * primi TMAX+1 valori di Dt
	 */
	for(t=0; t<TMAX+1; t++)
	{
		autocorrelation(corrDtstep, t, autocorr, N, STEPS);
		for(Dt=0; Dt<DMAX+1; Dt++)
		{
			fprintf(out_autocorr[Dt], "%d\t%e\n", t, autocorr[Dt]);
		}
	}

	for(i=0; i<DMAX+1; i++)
		fclose(out_autocorr[i]);
	for(i=0; i<N; i++)
	{
		free((clusterDt + i)->Vec);
		free((clusterSQDt + i)->Vec);
	}
	free(clusterDt);
	free(clusterSQDt);
	fclose(out_corr);
	fclose(out_action);
	fclose(out_deltaE);
	fclose(out_EM);
	fclose(out_x2gs);
	fclose(out_corrsq);
	fclose(out_deltaE_errors);
	fclose(out_EM_errors);

	exit(EXIT_SUCCESS);
}
