#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Compute the inverse chi-square CDF
int main (int argc, char **argv)
{
  if ( argc != 3 )
    {
      fprintf(stderr,"Incorrect arguments\n");
      return 1;
    }

  double x1 = atof(argv[1]);
  double x2 = atof(argv[2]);

  printf("%g\n",gsl_cdf_chisq_Pinv(x1,x2));

  return 0;
}
