#ifndef MATRIXH
#define MATRIXH

#include "stdio.h"
#include "complex.h"

typedef struct {
   Complex z11, z12, z21, z22;
} Matrix;

#ifdef USE_INLINE
#include "matrix.c"
#else
/* Matrix zero */
void mz(Matrix *Mptr);

/* Matrix identity */
void mi(Matrix *Mptr);

/* Matrix multiply */
Matrix mm(Matrix M1, Matrix M2);

/* Matrix distance */
double md(Matrix M1, Matrix M2);

/* Matrix inverse */
Matrix minv(Matrix M1);

/* Print matrix */
void pm(FILE *out, Matrix M);

/**
 * Matrix distance satisfying the triangle inequality
 * (albeit with two square roots in the computation!).
 */
double md_tri(Matrix M1, Matrix M2);

#endif

#endif
