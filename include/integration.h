#ifndef INTEGRATION_H
#define INTEGRATION_H

/* struttura integral:
 * contiene risultato (R) ed errore (E) di un integrale (valori float)
 */
typedef struct
{
	float R, E;
} integral;


/* struttura integral_dble:
 * contiene risultato (R) ed errore (E) di un integrale (valori double)
 */
typedef struct
{
	double R, E;
} integral_dble;


#ifndef INTEGRATE_C
extern double integration (double (*f)(double), double (*method)(double (*g)(double), double, double), double a, double b, int N);
#endif

#ifndef METHODS_C
extern float fsimpson (float (*f)(float), float a, float b);
extern double simpson (double (*f)(double), double a, double b);
extern float ftrapezoid (float (*f)(float), float a, float b);
extern double trapezoid (double (*f)(double), double a, double b);
extern double gaussianquad_5 (double (*f)(double), double a, double b);
#endif

#ifndef GAUSSMOM_C
extern void gaussian_moment (int m, void (*random_generator)(double *, int), double (*g)(double), int N, integral_dble * I);
extern void wgaussian_moment (int m, void (*wrandom_generator)(double *, int, double), double (*wg)(double, double), double A, int N, integral_dble * I);
#endif

#ifndef PIGRECO_C
extern double Pi_trivial(int N);
extern double Pi_buffon(int N, double L, double pigreco);
#endif

#endif
