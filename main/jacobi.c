
#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>
/*#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>*/
#include <gsl/gsl_eigen.h>
#include "extras.h"

int main (int argc, char *argv[])
{
	double data[] = {	0.0	,	1.0,
						-1.0,	0.0};
	
	gsl_matrix_view M = gsl_matrix_view_array (data, 2, 2);
	
	gsl_vector_complex *eigenvalues		= gsl_vector_complex_alloc (2);
	gsl_matrix_complex *eigenvectors	= gsl_matrix_complex_alloc (2,2);
	
	gsl_eigen_nonsymmv_workspace *W = gsl_eigen_nonsymmv_alloc (2);
	gsl_eigen_nonsymmv (&M.matrix, eigenvalues, eigenvectors, W);
	
	gsl_eigen_nonsymmv_free (W);
	
	int i;
	for (i = 0; i < 2; i++)
	{
		gsl_complex					eval_i = gsl_vector_complex_get (eigenvalues, i);
		gsl_vector_complex_view		evec_i = gsl_matrix_complex_column (eigenvectors, i);
		printf ("eigenvalue = %lf + i %lf\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));
		printf ("eigenvector = \n");
		gsl_vector_complex_fprintf (stdout,	&evec_i.vector, "%g");
	}

	gsl_vector_complex_free (eigenvalues);
	gsl_matrix_complex_free (eigenvectors);
	
	exit(EXIT_SUCCESS);
}
