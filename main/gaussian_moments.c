
/*******************************************************************************
 *
 *		File "gaussian_moments.c"
 *
 * Programma per il calcolo di momenti della gaussiana standard con metodi
 * Monte Carlo (integrazione con campionamento di importanza).
 *
 * Le funzioni utilizzate per il calcolo sono definite in:
 * _ ../modules/random (generatori di numeri casuali)
 *			gauss.c    symexpo.c   symrootexpo.c
 * _ ../modules/integration (integrazione Monte Carlo)
 *			gaussmom.c
 *
 ******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "random.h"
#include "integration.h"


int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		printf("\nParametri insufficienti. Impostare:\n");
		printf("\t1 _ min # di numeri random;\n");
		printf("\t2 _ max # di numeri random;\n");
		printf("\t3 _ fattore di scala logaritmica;\n");
		printf("\t4 _ massimo ordine di momento da calcolare.\n\n");
		
		exit(EXIT_FAILURE);
	}

	/* Eliminazione dati vecchi e creazione cartella per dati nuovi */
	system("rm -r MG");
	system("mkdir MG");
	
	int N_max, N_min, scale, N_points;
	int order, max_order;
	int choice;
	integral_dble *Moment;

	/* File di output (definito a seconda dei casi) */
	FILE *results;

	/* Stringhe per definire il nome del file di output */
	char *Path, *Type, *Order;
		Path = (char*)malloc(sizeof(char));
		Type = (char*)malloc(sizeof(char));
		Order = (char*)malloc(sizeof(char));
	const char *Ext = ".dat";

	N_min = atoi(argv[1]);
	N_max = atoi(argv[2]);
	scale = atoi(argv[3]);
	max_order = atoi(argv[4]);

	Moment = malloc(sizeof(integral_dble));
	
	/* Puntatori a funzioni utili:
	 * _ generator -> routine per generazione di numeri casuali secondo la
	 *		distribuzione specificata;
	 * _ distribution -> funzione corrispondente alla distribuzione dei numeri
	 *		prodotti da "generator".
	 */
	void (*generator)(double *, int);
	double (*distribution)(double);


	/*
	 * Calcolo e stampa dei momenti gaussiani di ordine da 1 a max_order
	 */
	for(order=1; order<max_order+1; order++)
	{
		for(choice=1; choice<4; choice++)
		{
			sprintf(Path, "MG/");
			sprintf(Order, "_%d", order);
			
			/*
			 * Scelta della distribuzione di numeri random da utilizzare
			 * e definizione del file di output
			 */
			switch (choice)
			{
				case 1:
					/* distribuzione gaussiana standardizzata */
					generator = stdgauss_dble;
					distribution = stdgaussdistr;
					sprintf(Type, "stdgauss");
					strcat(Path, Type);
					strcat(Path, Order);
					strcat(Path, Ext);
					results = fopen(Path,"w");
					break;

				case 2:
					/* distribuzione esponenziale */
					generator = symexpo_dble;
					distribution = symexpodistr;
					sprintf(Type, "expo");
					strcat(Path, Type);
					strcat(Path, Order);
					strcat(Path, Ext);
					results = fopen(Path,"w");
					break;

				case 3:
					/* distribuzione radice*esponenziale */
					generator = symrootexpo_dble;
					distribution = symrootexpodistr;
					sprintf(Type, "rootexpo");
					strcat(Path, Type);
					strcat(Path, Order);
					strcat(Path, Ext);
					results = fopen(Path,"w");
					break;
			}

			fprintf(results,"#\n# File contenente risultati per il momento gaussiano ");
			fprintf(results,"(gaussiana standardizzata) di ordine %d;\n", order);
			fprintf(results,"# risultati ottenuti con metodi Monte Carlo ");
			fprintf(results,"(campionamento di importanza con distribuzione \"%s\")\n", Type);
			fprintf(results,"# Le colonne sono:\n");
			fprintf(results,"# 1_# numeri random   2_Risultato   3_Errore (deviazione standard della media)\n#\n");
			
			printf("Momento di ordine %d; distribuzione \"%s\";\n", order, Type);
			fflush(stdout);

			/*
			 * Calcolo e stampa del momento gaussiano
			 * (ciclo sul numero di punti random da generare)
			 */
			for(N_points=N_min; N_points<N_max+1; N_points *= scale)
			{
				printf("\t%d\tpunti...\n",N_points);
				gaussian_moment(order, generator, distribution, N_points, Moment);
				fprintf(results, "%d\t%e\t%e\n", N_points, (*Moment).R, (*Moment).E);
			}
			printf("\t\tCompletata!\n");
		}
	}

	fclose(results);
	printf("\n\nRisultati salvati in \"MG\";\n");
	printf("i files sono \"<distribuzione usata>_<ordine del momento>.dat\"\n\n");
	sleep(2);
	
	
	exit(EXIT_SUCCESS);
}
