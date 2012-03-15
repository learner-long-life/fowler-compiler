#ifndef COMPLEXH
#define COMPLEXH

typedef struct {
   double x, y;
} Complex;

#ifdef USE_INLINE
#include "complex.c"
#else
/* Complex multiply */
Complex cm(Complex z1, Complex z2);

/* Complex add */
Complex ca(Complex z1, Complex z2);

/* Complex subtract */
Complex cs(Complex z1, Complex z2);

/* Complex conjugate */
Complex cc(Complex z);

/* Complex absolute value */
double my_cabs(Complex z);

/* Complex zero test */
int cz(Complex z);

#endif

#endif
