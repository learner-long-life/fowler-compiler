#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include "complex.h"
#include "matrix.h"

/* Number gates */
#define NG 26

/* Gate numbers */
#define I    0
#define H    1
#define X    2
#define Z    3
#define S    4
#define Sd   5
#define HX   6
#define HZ   7
#define HS   8
#define HSd  9
#define XZ   10
#define XS   11
#define XSd  12
#define SH   13
#define SdH  14
#define HXZ  15
#define HXS  16
#define HXSd 17
#define HSH  18
#define HSdH 19
#define XSH  20
#define XSdH 21
#define SHSd 22
#define SdHS 23
#define T    24
#define Td   25

/* Floating point format */
#define FRAC_BITS 35
#define TOTAL_BITS 37

double sqrt2o2;
Matrix gate_list[NG];

void init_matrices() {
    sqrt2o2 = sqrt(2)/2;

    /* Identity gate */
    mi(&gate_list[I]);

    /* Hadamard gate */
    mz(&gate_list[H]);
    gate_list[H].z11.x = sqrt2o2;
    gate_list[H].z12.x = sqrt2o2;
    gate_list[H].z21.x = sqrt2o2;
    gate_list[H].z22.x = -sqrt2o2;

    /* X gate */
    mz(&gate_list[X]);
    gate_list[X].z12.x = 1;
    gate_list[X].z21.x = 1;

    /* Z gate */
    mz(&gate_list[Z]);
    gate_list[Z].z11.x = 1;
    gate_list[Z].z22.x = -1;

    /* S gate */
    mz(&gate_list[S]);
    gate_list[S].z11.x = 1;
    gate_list[S].z22.y = 1;

    /* Sd gate */
    mz(&gate_list[Sd]);
    gate_list[Sd].z11.x = 1;
    gate_list[Sd].z22.y = -1;

    /* HX gate */
    mz(&gate_list[HX]);
    gate_list[HX] = mm(gate_list[X],gate_list[H]);

    /* HZ gate */
    mz(&gate_list[HZ]);
    gate_list[HZ] = mm(gate_list[Z],gate_list[H]);

    /* HS gate */
    mz(&gate_list[HS]);
    gate_list[HS] = mm(gate_list[S],gate_list[H]);

    /* HSd gate */
    mz(&gate_list[HSd]);
    gate_list[HSd] = mm(gate_list[Sd],gate_list[H]);

    /* XZ gate */
    mz(&gate_list[XZ]);
    gate_list[XZ] = mm(gate_list[Z],gate_list[X]);

    /* XS gate */
    mz(&gate_list[XS]);
    gate_list[XS] = mm(gate_list[S],gate_list[X]);

    /* XSd gate */
    mz(&gate_list[XSd]);
    gate_list[XSd] = mm(gate_list[Sd],gate_list[X]);

    /* SH gate */
    mz(&gate_list[SH]);
    gate_list[SH] = mm(gate_list[H],gate_list[S]);

    /* SdH gate */
    mz(&gate_list[SdH]);
    gate_list[SdH] = mm(gate_list[H],gate_list[Sd]);

    /* HXZ gate */
    mz(&gate_list[HXZ]);
    gate_list[HXZ] = mm(gate_list[Z],mm(gate_list[X],gate_list[H]));

    /* HXS gate */
    mz(&gate_list[HXS]);
    gate_list[HXS] = mm(gate_list[S],mm(gate_list[X],gate_list[H]));

    /* HXSd gate */
    mz(&gate_list[HXSd]);
    gate_list[HXSd] = mm(gate_list[Sd],mm(gate_list[X],gate_list[H]));

    /* HSH gate */
    mz(&gate_list[HSH]);
    gate_list[HSH] = mm(gate_list[H],mm(gate_list[S],gate_list[H]));

    /* HSdH gate */
    mz(&gate_list[HSdH]);
    gate_list[HSdH] = mm(gate_list[H],mm(gate_list[Sd],gate_list[H]));

    /* XSH gate */
    mz(&gate_list[XSH]);
    gate_list[XSH] = mm(gate_list[H],mm(gate_list[S],gate_list[X]));

    /* XSdH gate */
    mz(&gate_list[XSdH]);
    gate_list[XSdH] = mm(gate_list[H],mm(gate_list[Sd],gate_list[X]));

    /* SHSd gate */
    mz(&gate_list[SHSd]);
    gate_list[SHSd] = mm(gate_list[Sd],mm(gate_list[H],gate_list[S]));

    /* SdHS gate */
    mz(&gate_list[SdHS]);
    gate_list[SdHS] = mm(gate_list[S],mm(gate_list[H],gate_list[Sd]));

    /* T gate */
    mz(&gate_list[T]);
    gate_list[T].z11.x = 1;
    gate_list[T].z22.x = sqrt2o2;
    gate_list[T].z22.y = sqrt2o2;

    /* Td gate */
    mz(&gate_list[Td]);
    gate_list[Td].z11.x = 1;
    gate_list[Td].z22.x = sqrt2o2;
    gate_list[Td].z22.y = -sqrt2o2;
}

double to_fixed_point(double number) {
    return (number * ((int64_t)1 << FRAC_BITS));
}

FILE *fp;
int index;

void print_single_entry(double entry) {
  fprintf(fp, "\t%d\t:\t%.0f;\n", index, to_fixed_point(entry));
  index++;
}

void print_complex_number(Complex *num) {
  print_single_entry(num->x);
  print_single_entry(num->y);
}

void print_matrix(Matrix *mtx) {
  print_complex_number(&(mtx->z11));
  print_complex_number(&(mtx->z12));
  print_complex_number(&(mtx->z21));
  print_complex_number(&(mtx->z22));
}

int print_matrices() {
    fp = fopen("gates_37bit.mif", "wb");
    if (fp == NULL) {
      fprintf(stderr, "Failed to open output file.\n");
      return 1;
    }
    index = 0;

    fprintf(fp, "WIDTH=%d;\n", TOTAL_BITS);
    fprintf(fp, "DEPTH=%d;\n", (NG - 1) * 8);
    fprintf(fp, "ADDRESS_RADIX=DEC;\n");
    fprintf(fp, "DATA_RADIX=DEC;\n");
    fprintf(fp, "CONTENT BEGIN\n");

    int i;
    for (i = 1; i < NG; i++) { print_matrix(gate_list + i); }

    fprintf(fp, "END\n");

    fclose(fp);
    return 0;
}

int main() {
    init_matrices();
    return print_matrices();
}

