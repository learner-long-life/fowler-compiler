#ifndef MATRIXH
#define MATRIXH

#include "stdio.h"
#include "complex.h"

typedef struct {
   Complex z11, z12, z21, z22;
} Matrix;

/* Matrix zero */
void mz(Matrix *Mptr);

/* Matrix identity */
void mi(Matrix *Mptr);

/* Matrix multiply */
Matrix mm(Matrix M1, Matrix M2);

/* Matrix distance */
double md(Matrix M1, Matrix M2);

/* Print matrix */
void pm(FILE *out, Matrix M);

#endif
