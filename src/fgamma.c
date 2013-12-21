#include <Rmath.h>      /* to define fgamma, etc */
#include <R_ext/RS.h>   /* to define F77_NAME */

double F77_NAME(fgamma)(double *x)
{
  return gammafn(*x);
}
