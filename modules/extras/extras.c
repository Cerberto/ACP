
/*******************************************************************************
 *
 *		File "extras.c"
 *
 * File contenente routines extra:
 * _ factorial -> calcolo del fattoriale.
 *
 ******************************************************************************/

#define EXTRAS_C

#include <stdio.h>

int factorial(int n)
{
	int res = 1;
	int i;
	for(i=1; i<n+1; i++)
		res *= i;
	return res;
}

void pippo()
{
	printf("\tPippo!\n");
	fflush(stdout);
}
