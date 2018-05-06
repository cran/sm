#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(diag_cov_bin_fun)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(full_cov_bin_fun)(void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"diag_cov_bin_fun", (DL_FUNC) &F77_NAME(diag_cov_bin_fun), 7},
    {"full_cov_bin_fun", (DL_FUNC) &F77_NAME(full_cov_bin_fun), 7},
    {NULL, NULL, 0}
};

void R_init_sm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
