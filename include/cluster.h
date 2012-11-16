
#ifndef CLUSTER_H
#define CLUSTER_H


typedef struct
{
	double *Vec;	/* Array of data */
	int Dim;		/* Array length */
	double Mean;	/* Mean value on the array */
	double Sigma;	/* Variance of the mean */
} cluster;


#ifndef JACKKNIFE_C
extern void clusterJK(cluster *C);
extern void cluster_init(cluster *C, int dim);
extern cluster functionJK (double (*f)(double), cluster *X);
#endif

#endif

