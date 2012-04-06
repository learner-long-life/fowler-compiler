#ifndef MATRIXH
#define MATRIXH

#include "stdio.h"

#ifndef PAULI_BASIS
#include "complex.h"
#endif

#ifdef PAULI_BASIS
class Matrix {
public:
  double a[4];
  Matrix();
  Matrix(double a0, double a1, double a2, double a3);
  double &operator[](const int index);
  double operator[](const int index) const;
  double dot(Matrix &b) const;
  Matrix operator*(const double n) const;
  Matrix operator*(const Matrix m) const;
};
#else
typedef struct {
   Complex z11, z12, z21, z22;
} Matrix;
#endif

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
