#include <math.h>
#include "complex.h"

#ifdef COMPLEX_SIMD
typedef double v2df __attribute__ ((vector_size(16))); // vector of 2 double floats

union f2vec {
  double d[2];
  v2df v;
};
#endif

#ifdef USE_INLINE
inline
#endif
Complex cm(Complex z1, Complex z2) {
   Complex z;

#ifdef COMPLEX_SIMD
   union f2vec a, b, c;
   a.d[0] = z1.x;
   a.d[1] = z1.y;
   b.d[0] = z2.x;
   b.d[1] = z2.y;
   c.v = a.v * b.v;
   z.x = c.d[0] - c.d[1];

   b.d[0] = z2.y;
   b.d[1] = z2.x;
   c.v = a.v * b.v;
   z.y = c.d[0] + c.d[1];
#else
   z.x = z1.x * z2.x - z1.y * z2.y;
   z.y = z1.x * z2.y + z1.y * z2.x;
#endif

   return z;
}

#ifdef USE_INLINE
inline
#endif
Complex ca(Complex z1, Complex z2) {
   Complex z;

   z.x = z1.x + z2.x;
   z.y = z1.y + z2.y;

   return z;
}

#ifdef USE_INLINE
inline
#endif
Complex cs(Complex z1, Complex z2) {
   Complex z;

   z.x = z1.x - z2.x;
   z.y = z1.y - z2.y;

   return z;
}

#ifdef USE_INLINE
inline
#endif
Complex cc(Complex z) {
   z.y = -z.y;

   return z;
}

#ifdef USE_INLINE
inline
#endif
double my_cabs(Complex z) {
   return z.x * z.x + z.y * z.y;
}

#ifdef USE_INLINE
inline
#endif
int cz(Complex z) {
   if (z.x==0 && z.y==0) return 0;
   return 0;
}

