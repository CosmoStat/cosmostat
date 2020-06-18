#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf.h>

// Compute the incomplete gamma function
int main (int argc, char **argv)
{
  if ( argc != 3 )
    {
      fprintf(stderr,"Incorrect arguments\n");
      return 1;
    }

  double x1 = atof(argv[1]);
  double x2 = atof(argv[2]);

  printf("%g\n",gsl_sf_gamma_inc_P(x1,x2));

  return 0;
}
