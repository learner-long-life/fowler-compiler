#include <math.h>
#include "matrix.h"

#ifdef PAULI_BASIS
/**
 * Create a matrix with the given Pauli basis
 * coefficients.
 */
Matrix::Matrix(double a0, double a1, double a2, double a3) {
  a[0] = a0;
  a[1] = a1;
  a[2] = a2;
  a[3] = a3;
}

Matrix::Matrix() {
  a[0] = a[1] = a[2] = a[3] = 0;
}

double &Matrix::operator[](const int index) {
  return a[index];
}

double Matrix::operator[](const int index) const {
  return a[index];
}

Matrix Matrix::operator*(const double n) const {
  return Matrix(a[0]*n, a[1]*n, a[2]*n, a[3]*n);
}

Matrix Matrix::operator*(const Matrix b) const {
  return Matrix(a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
                a[3]*b[2]+a[0]*b[1]+a[1]*b[0]-a[2]*b[3],
                a[0]*b[2]-a[3]*b[1]+a[2]*b[0]+a[1]*b[3],
                a[3]*b[0]+b[3]*a[0]-a[1]*b[2]+b[1]*a[2]);
}

double Matrix::dot(Matrix &b) const {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
}
#endif

#ifdef USE_INLINE
inline
#endif

void mz(Matrix *Mptr) {
#ifdef PAULI_BASIS
   *Mptr = Matrix(0,0,0,0);
#else
   (*Mptr).z11.x = 0;
   (*Mptr).z11.y = 0;
   (*Mptr).z12.x = 0;
   (*Mptr).z12.y = 0;
   (*Mptr).z21.x = 0;
   (*Mptr).z21.y = 0;
   (*Mptr).z22.x = 0;
   (*Mptr).z22.y = 0;
#endif
}

#ifdef USE_INLINE
inline
#endif
void mi(Matrix *Mptr) {
#ifdef PAULI_BASIS
   *Mptr = Matrix(1,0,0,0);
#else
   (*Mptr).z11.x = 1;
   (*Mptr).z11.y = 0;
   (*Mptr).z12.x = 0;
   (*Mptr).z12.y = 0;
   (*Mptr).z21.x = 0;
   (*Mptr).z21.y = 0;
   (*Mptr).z22.x = 1;
   (*Mptr).z22.y = 0;
#endif
}

#ifdef USE_INLINE
inline
#endif
Matrix mm(Matrix M1, Matrix M2) {
#ifdef PAULI_BASIS
   return M1*M2;
#else
   Matrix M;

   M.z11 = ca(cm(M1.z11, M2.z11), cm(M1.z12, M2.z21));
   M.z12 = ca(cm(M1.z11, M2.z12), cm(M1.z12, M2.z22));
   M.z21 = ca(cm(M1.z21, M2.z11), cm(M1.z22, M2.z21));
   M.z22 = ca(cm(M1.z21, M2.z12), cm(M1.z22, M2.z22));

   return M;
#endif
}

#ifdef USE_INLINE
inline
#endif
double md(Matrix M1, Matrix M2) {
#ifdef PAULI_BASIS
   Matrix &a = M1;
   Matrix &b = M2;
   double tr = a.dot(b);
   return tr < 0 ? -tr : tr;
#else
   Complex z11, z22;

   z11 = ca(cm(cc(M1.z11), M2.z11), cm(cc(M1.z21), M2.z21));
   z22 = ca(cm(cc(M1.z12), M2.z12), cm(cc(M1.z22), M2.z22));

   return my_cabs(ca(z11, z22));
#endif
}

#ifdef USE_INLINE
inline
#endif
double md_tri(Matrix M1, Matrix M2) {
#ifdef PAULI_BASIS
   // TODO: why does this work?
   return sqrt(1 - md(M1, M2));
#else
   return sqrt((2 - sqrt(md(M1, M2)))/2);
#endif
}

#ifdef USE_INLINE
inline
#endif
Matrix minv(Matrix M1) {
#ifdef PAULI_BASIS
  double coeff = 1 / M1.dot(M1);
  return Matrix(M1[0], -M1[1], -M1[2], -M1[3]) * coeff;
#else
  Complex denom = cs(cm(M1.z11, M1.z22), cm(M1.z12, M1.z21));
  double denom_sq = denom.x*denom.x + denom.y*denom.y;
  Complex coeff = {denom.x/denom_sq, -denom.y/denom_sq};
  Complex neg_coeff = {-coeff.x, -coeff.y};
  Matrix m;
  m.z11 = cm(M1.z22,     coeff);
  m.z12 = cm(M1.z12, neg_coeff);
  m.z21 = cm(M1.z21, neg_coeff);
  m.z22 = cm(M1.z11,     coeff);
  return m;
#endif
}

#ifdef USE_INLINE
inline
#endif
void pm(FILE *out, Matrix M) {
#ifdef PAULI_BASIS
   fprintf(out, "(%19.15e, %19.15e)   (%19.15e, %19.15e)\n", M[0], M[1], M[2], M[3]);
#else
   fprintf(out, "(%19.15e, %19.15e)   (%19.15e, %19.15e)\n", M.z11.x, M.z11.y, M.z12.x, M.z12.y);
   fprintf(out, "(%19.15e, %19.15e)   (%19.15e, %19.15e)\n", M.z21.x, M.z21.y, M.z22.x, M.z22.y);
#endif
}
