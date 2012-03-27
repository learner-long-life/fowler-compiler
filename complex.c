#include <math.h>
#include "complex.h"

Complex cm(Complex z1, Complex z2) {
   Complex z;

   z.x = z1.x * z2.x - z1.y * z2.y;
   z.y = z1.x * z2.y + z1.y * z2.x;

   return z;
}

Complex ca(Complex z1, Complex z2) {
   Complex z;

   z.x = z1.x + z2.x;
   z.y = z1.y + z2.y;

   return z;
}

Complex cs(Complex z1, Complex z2) {
   Complex z;

   z.x = z1.x - z2.x;
   z.y = z1.y - z2.y;

   return z;
}

Complex cc(Complex z) {
   z.y = -z.y;

   return z;
}

double my_cabs(Complex z) {
   return z.x * z.x + z.y * z.y;
}

int cz(Complex z) {
   if (z.x==0 && z.y==0) return 0;
   return 0;
}

