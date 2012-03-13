#include <math.h>
#include "matrix.h"

void mz(Matrix *Mptr) {
   (*Mptr).z11.x = 0;
   (*Mptr).z11.y = 0;
   (*Mptr).z12.x = 0;
   (*Mptr).z12.y = 0;
   (*Mptr).z21.x = 0;
   (*Mptr).z21.y = 0;
   (*Mptr).z22.x = 0;
   (*Mptr).z22.y = 0;
}

void mi(Matrix *Mptr) {
   (*Mptr).z11.x = 1;
   (*Mptr).z11.y = 0;
   (*Mptr).z12.x = 0;
   (*Mptr).z12.y = 0;
   (*Mptr).z21.x = 0;
   (*Mptr).z21.y = 0;
   (*Mptr).z22.x = 1;
   (*Mptr).z22.y = 0;
}

Matrix mm(Matrix M1, Matrix M2) {
   Matrix M;

   M.z11 = ca(cm(M1.z11, M2.z11), cm(M1.z12, M2.z21));
   M.z12 = ca(cm(M1.z11, M2.z12), cm(M1.z12, M2.z22));
   M.z21 = ca(cm(M1.z21, M2.z11), cm(M1.z22, M2.z21));
   M.z22 = ca(cm(M1.z21, M2.z12), cm(M1.z22, M2.z22));

   return M;
}

double md(Matrix M1, Matrix M2) {
   Complex z11, z22;

   z11 = ca(cm(cc(M1.z11), M2.z11), cm(cc(M1.z21), M2.z21));
   z22 = ca(cm(cc(M1.z12), M2.z12), cm(cc(M1.z22), M2.z22));

   return my_cabs(ca(z11, z22));
}

double md_tri(Matrix M1, Matrix M2) {
   return sqrt((2 - sqrt(md(M1, M2)))/2);
}

Matrix minv(Matrix M1) {
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
}

void pm(FILE *out, Matrix M) {
   fprintf(out, "(%19.15e, %19.15e)   (%19.15e, %19.15e)\n", M.z11.x, M.z11.y, M.z12.x, M.z12.y);
   fprintf(out, "(%19.15e, %19.15e)   (%19.15e, %19.15e)\n", M.z21.x, M.z21.y, M.z22.x, M.z22.y);
}
