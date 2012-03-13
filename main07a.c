#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "bin.h"
#include "complex.h"
#include "matrix.h"

#ifdef BENCHMARK
#include "time.h"
#endif

/* Number gates */
#define NG 26

/* Array sizes */
#define SMALL 100
#define MEDIUM 10000
#define BIG 10000000

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

void print_gate(FILE *out, int gate);
void print_product(FILE *out, int *product, int most_significant);
void print_unique_matrices(FILE *out);
void print_unique_product_lists(FILE *out);
void print_product_check_tree(FILE *out);

int product_less_than(int *product1, int *product2, int most_significant);
int unique_greater_than_product(int unique);
int increment_product(int start);
int test_product(int start, int end);
int overwrite_product(int start, int list);
int skip_product(int last_most_significant); 
int skip_product_print(int last_most_significant); 
int skip_product_long(int last_most_significant, int width); 
int is_unique_product();
void calculate_product(int last_most_significant);
void add_matrix_to_unique();
void insert_product_in_list();
void insert_product_in_tree();
void add_matrix_to_structures();

double Pi, sqrt2o2, epsilon;
int *product;
int *unique_product_lists;
Matrix *unique_matrices;
int *product_check_tree;
int most_significant, free_list, num_unique, free_node;
Matrix U, U1, U2, U3, U4, gate_list[NG];

FILE *out;

#ifdef BIN
BinSet seq_bins(0, 1, 10000);
double gate_accuracy = .112;

void print_sequence(FILE* out, int index) {
   int end = index + unique_product_lists[index];
   for (int i = index + 1; i <= end; i++) {
     print_gate(out, unique_product_lists[i]);
   }
   fprintf(out, "\n");
}
#endif

//#ifdef DISTANCES
//double *distances;
//int next_dist;
//#endif

int main() {
   int last_most_significant, input_format, numerator, denominator;
   int n, width;
   double dist, temp_dist, x, y;
   Matrix G;
   FILE *in = (FILE *)fopen("in", "r");
   out = (FILE *)fopen("out", "w");

   /* Initialise global variables */
   Pi = acos(-1);
   sqrt2o2 = sqrt(2)/2;
   epsilon=1e-10;

   /* number of unique products to find */
   fscanf(in, "%d", &n);

   /* width of product check tree */
   fscanf(in, "%d", &width);
   fprintf(out, "number of unique products to find n: %d, width: %d\n", n, width);

   //#ifdef DISTANCES
   //distances = (double*)calloc(n, sizeof(double));
   //next_dist = 0;
   //#endif

   product=(int *)calloc(SMALL, sizeof(int));
   if (product==NULL) exit(0);
   most_significant=0;
   fprintf(out, "Size of product = %d\n", SMALL);

   /* unique products are generally short, so 25*n is generous memory */
   unique_product_lists=(int *)calloc(25*n, sizeof(int));
   if (unique_product_lists==NULL) exit(0);
   free_list=2;
   fprintf(out, "Size of unique_product_lists = %d\n", 25*n);

   unique_matrices=(Matrix *)calloc(n, sizeof(Matrix));
   if (unique_matrices==NULL) exit(0);
   mi(&unique_matrices[0]);
   num_unique=1;
   fprintf(out, "Size of unique_matrices = %d\n", n);

   product_check_tree=(int *)calloc((NG+2)*n, sizeof(int));
   if (product_check_tree==NULL) exit(0);
   free_node=NG+1;
   fprintf(out, "Size of product_check_tree = %d\n", (NG+2)*n);

   mi(&U);
   mi(&U1);
   mi(&U2);
   mi(&U3);
   mi(&U4);

   /* Gate to approximate */
   mz(&G);
   fscanf(in, "%d", &input_format);
   fprintf(out, "input_format: %d\n", input_format);
   if (input_format==0) {
      fscanf(in, "%lf", &G.z11.x);
      fscanf(in, "%lf", &G.z11.y);
      fscanf(in, "%lf", &G.z12.x);
      fscanf(in, "%lf", &G.z12.y);
      fscanf(in, "%lf", &G.z21.x);
      fscanf(in, "%lf", &G.z21.y);
      fscanf(in, "%lf", &G.z22.x);
      fscanf(in, "%lf", &G.z22.y);
   }
   else if (input_format==1) {
      fscanf(in, "%d", &numerator);
      fscanf(in, "%d", &denominator);
      fprintf(out, "numerator: %d, denominator: %d\n", numerator, denominator);
      G.z11.x=1;
      x = Pi*numerator/denominator;
      G.z22.x=cos(x);
      G.z22.y=sin(x);
   }
   fclose(in);

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

#ifdef BIN
   // use the triangle-inequality distance measure
   dist = md_tri(gate_list[I], G);
#else
   dist = md(gate_list[I],G);
#endif
   fprintf(out, "dist=%17.13e\n", dist);

   fprintf(stderr, "G-1 * G should be the identity:\n");
   pm(stderr, mm(minv(gate_list[Td]), gate_list[Td]));

#ifdef BENCHMARK
   printf("Beginning enumeration...\n", num_unique);
   struct timespec start, end;
   clock_gettime(CLOCK_MONOTONIC, &start);
   int last_ms = 1;
#endif

   while (most_significant+1<=width && num_unique < n) {
      last_most_significant=increment_product(0);
      last_most_significant = skip_product(last_most_significant);
      calculate_product(last_most_significant);

      if (is_unique_product()) {
#ifdef BIN
         temp_dist = md_tri(U1,G);
         int seq_index = seq_bins.contains(U1, temp_dist, dist - epsilon, true,
                                           dist);
         if (seq_index != -1) {
            fprintf(out, "Found 'meet in the middle' sequence with distance %.10f: \n",
                    dist);
            print_product(out, product, most_significant);
            print_sequence(out, seq_index);
            fflush(out);
         }
         // TODO: replicate seq_bins.contains in the second stage check.
         Matrix inv = mm(G, minv(U1));
         seq_bins.insert(md_tri(inv, G), inv, free_list, most_significant + 1, true);
#else
         temp_dist=md(U1,G);
#endif
         add_matrix_to_structures();
         #ifdef DISTANCES
         printf("Dist: %.9f\n", temp_dist);
         #endif
#ifdef BIN
         // the distance decreases, so signs get swapped
         if (temp_dist < dist + epsilon) {
            if (temp_dist < dist - epsilon) {
#else
         if (temp_dist>dist-epsilon) {
            if (temp_dist>dist+epsilon) {
#endif
               dist=temp_dist;
               fprintf(out, "dist=%17.13e\n", dist);
            }
            print_product(out, product, most_significant);
            fflush(out);
         }
      }
#ifdef BENCHMARK
      if (last_ms != most_significant) {
         last_ms = most_significant;
         clock_gettime(CLOCK_MONOTONIC, &end);
         double start_ns = start.tv_sec * 1000000000. + start.tv_nsec;
         double end_ns   =   end.tv_sec * 1000000000. +   end.tv_nsec;
         double duration = (end_ns-start_ns)/1000000000.;
         printf("Time to calculate %d sequences to length %d: %.9f\n",
                num_unique, last_ms, duration);
      }
#endif
   }

   fprintf(out, "num_unique: %d\n", num_unique);
#ifdef BENCHMARK
   {
      clock_gettime(CLOCK_MONOTONIC, &end);
      double start_ns = start.tv_sec * 1000000000. + start.tv_nsec;
      double end_ns   =   end.tv_sec * 1000000000. +   end.tv_nsec;
      double duration = (end_ns-start_ns)/1000000000.;
      printf("num_unique: %d\n", num_unique);
      printf("Elapsed time: %.9f\n", duration);
   }
#endif
#ifdef BIN
   seq_bins.delete_short_sequences(width);
#endif
#ifdef DISTANCES
   return 0;
#endif

   /*
   print_unique_product_lists(out);
   print_product_check_tree(out);
   */

   while (1) {
      last_most_significant=increment_product(0);
      last_most_significant = skip_product_long(last_most_significant, width);
      calculate_product(last_most_significant);
#ifdef BIN
      temp_dist = md_tri(U1,G);
      int seq_index = seq_bins.contains(U1, temp_dist, dist - epsilon, true,
                                        dist);
      if (seq_index != -1) {
         fprintf(out, "Found 'meet in the middle' sequence with distance %.10f: \n",
                 dist);
         print_product(out, product, most_significant);
         print_sequence(out, seq_index);
         fflush(out);
      }

      // the distance decreases, so signs get swapped
      if (temp_dist < dist + epsilon) {
         if (temp_dist < dist - epsilon) {
#else
      temp_dist=md(U1,G);
      if (temp_dist>dist-epsilon) {
         if (temp_dist>dist+epsilon) {
#endif
            dist=temp_dist;
            fprintf(out, "dist=%17.13e\n", dist);
         }
         print_product(out, product, most_significant);
         fflush(out);
      }
#ifdef BENCHMARK
      if (last_ms != most_significant) {
         last_ms = most_significant;
         clock_gettime(CLOCK_MONOTONIC, &end);
         double start_ns = start.tv_sec * 1000000000. + start.tv_nsec;
         double end_ns   =   end.tv_sec * 1000000000. +   end.tv_nsec;
         double duration = (end_ns-start_ns)/1000000000.;
         printf("Time to calculate %d sequences to length %d: %.9f\n",
                num_unique, last_ms, duration);
      }
#endif
   }

   fclose(out);
   free(product);
   free(unique_product_lists);
   free(unique_matrices);
   free(product_check_tree);

   return 0;
}

void print_gate(FILE *out, int gate) {
   if (gate==I) fprintf(out, "I");
   if (gate==H) fprintf(out, "H");
   if (gate==X) fprintf(out, "X");
   if (gate==Z) fprintf(out, "Z");
   if (gate==S) fprintf(out, "S");
   if (gate==Sd) fprintf(out, "Sd");
   if (gate==HX) fprintf(out, "(HX)");
   if (gate==HZ) fprintf(out, "(HZ)");
   if (gate==HS) fprintf(out, "(HS)");
   if (gate==HSd) fprintf(out, "(HSd)");
   if (gate==XZ) fprintf(out, "(XZ)");
   if (gate==XS) fprintf(out, "(XS)");
   if (gate==XSd) fprintf(out, "(XSd)");
   if (gate==SH) fprintf(out, "(SH)");
   if (gate==SdH) fprintf(out, "(SdH)");
   if (gate==HXZ) fprintf(out, "(HXZ)");
   if (gate==HXS) fprintf(out, "(HXS)");
   if (gate==HXSd) fprintf(out, "(HXSd)");
   if (gate==HSH) fprintf(out, "(HSH)");
   if (gate==HSdH) fprintf(out, "(HSdH)");
   if (gate==XSH) fprintf(out, "(XSH)");
   if (gate==XSdH) fprintf(out, "(XSdH)");
   if (gate==SHSd) fprintf(out, "(SHSd)");
   if (gate==SdHS) fprintf(out, "(SdHS)");
   if (gate==T) fprintf(out, "T");
   if (gate==Td) fprintf(out, "Td");
}

void print_product(FILE *out, int *product, int most_significant) {
   int i;

   for (i=most_significant;i>=0;i--) print_gate(out, product[i]);
   fprintf(out, "\n");
}

void print_unique_matrices(FILE *out) {
   int i;

   for (i=0; i<num_unique; i++) {
      pm(out, unique_matrices[i]);
      fprintf(out, "\n");
   }
}

void print_unique_product_lists(FILE *out) {
   int i, j;

   fprintf(out, "0, 1, I\n");
   i=2;
   while (i<free_list) {
      fprintf(out, "%d, %d, ", i, unique_product_lists[i]);
      for (j=1; j<=unique_product_lists[i]; j++) print_gate(out, unique_product_lists[i+j]);
      fprintf(out, "\n");
      i+=unique_product_lists[i]+1;
   }
}

void print_product_check_tree(FILE *out) {
   int i, j;

   i=0;
   while (i<free_node) {
      fprintf(out, "%4d", i);
      for (j=0; j<NG+1; j++) fprintf(out, "%4d ", product_check_tree[i+j]);
      fprintf(out, "\n");
      i+=NG+1;
   }
}

int product_less_than(int *product1, int *product2, int most_significant) {
   int i;

   for (i=most_significant; i>=0; i--) {
      if (product1[i]<product2[i]) return 1;
      if (product1[i]>product2[i]) return 0;
   }

   return 0;
}

int unique_greater_than_product(int unique) {
   int i;

   for (i=0; i<=most_significant; i++) {
      if (unique_product_lists[unique+i]<product[most_significant-i]) return 0;
      if (unique_product_lists[unique+i]>product[most_significant-i]) return 1;
   }

   return 1;
}

int increment_product(int start) {
   int i=start;
   
   for (i=0; i<start; i++) product[i]=1;
   
   i=start;
   product[i]++;
   while (product[i]>NG-1) {
      product[i]=1;
      i++;
      product[i]++;
   }

   if (i>most_significant) {
      most_significant=i;
      fprintf(out, "gate=%d\n", most_significant+1);
      fflush(out);
   }

   return i;
}

/* start >= end */
int test_product(int start, int end) {
   int i, current_position;

   current_position=0;
   for (i=start; i>=end; i--) {
      current_position+=product[i];
      current_position=product_check_tree[current_position];
      if (current_position<0) {
         if (unique_product_lists[-current_position]>start-i+1) return 1;
         else return current_position;
      }
   }

   return 0;
}

int overwrite_product(int start, int list) {
   int i, length, temp_lms;

   length=unique_product_lists[list];
   temp_lms=start-length+1;
   for (i=0; i<temp_lms; i++) product[i]=1;
   for (i=0; i<length; i++) {
      if (product[start-length+1+i]!=unique_product_lists[list+length-i])
         temp_lms=start-length+1+i;
      product[start-length+1+i]=unique_product_lists[list+length-i];
   }

   return temp_lms;
}

/* Returns new value of last_most_significant */
int skip_product(int last_most_significant) {
   int test, temp_lms;

   if (last_most_significant!=0) goto test_upper;

   test_lower:
   test=test_product(most_significant-1, 0);
   if (test==1) {
      temp_lms=increment_product(most_significant);
      if (temp_lms>last_most_significant) last_most_significant=temp_lms;
      goto test_upper;
   }
   if (test==0) return last_most_significant;
   temp_lms=overwrite_product(most_significant-1, -test);
   if (temp_lms>last_most_significant) last_most_significant=temp_lms;
   if (temp_lms==0) return last_most_significant;

   test_upper:
   test=test_product(most_significant, 1);
   if (test==1) {
      temp_lms=increment_product(most_significant+1);
      if (temp_lms>last_most_significant) last_most_significant=temp_lms;
      goto test_upper;
   }
   if (test==0) goto test_lower;
   temp_lms=overwrite_product(most_significant, -test);
   if (temp_lms>last_most_significant) last_most_significant=temp_lms;
   goto test_upper;

   return 0;
}

/* Returns new value of last_most_significant */
int skip_product_print(int last_most_significant) {
   int test, temp_lms;

   if (last_most_significant!=0) goto test_upper;

   test_lower:
   test=test_product(most_significant-1, 0);
   fprintf(out, "After testing lower, test=%d\n", test);
         fflush(out);
   if (test==1) {
      temp_lms=increment_product(most_significant);
      if (temp_lms>last_most_significant) last_most_significant=temp_lms;
      fprintf(out, "After major increment, lms=%d\n", last_most_significant);
      print_product(out, product, most_significant);
         fflush(out);
      goto test_upper;
   }
   if (test==0) return last_most_significant;
   temp_lms=overwrite_product(most_significant-1, -test);
   fprintf(out, "After overwriting product, temp_lms=%d\n", temp_lms);
   print_product(out, product, most_significant);
         fflush(out);
   if (temp_lms>last_most_significant) last_most_significant=temp_lms;
   if (temp_lms==0) return last_most_significant;

   test_upper:
   test=test_product(most_significant, 1);
   fprintf(out, "After testing upper, test=%d\n", test);
         fflush(out);
   if (test==1) {
      temp_lms=increment_product(most_significant+1);
      if (temp_lms>last_most_significant) last_most_significant=temp_lms;
      fprintf(out, "After major increment, lms=%d\n", last_most_significant);
      print_product(out, product, most_significant);
         fflush(out);
      goto test_upper;
   }
   if (test==0) goto test_lower;
   temp_lms=overwrite_product(most_significant, -test);
   fprintf(out, "After overwriting product, temp_lms=%d\n", temp_lms);
   print_product(out, product, most_significant);
         fflush(out);
   if (temp_lms>last_most_significant) last_most_significant=temp_lms;
   goto test_upper;

   return 0;
}

/* Returns new value of last_most_significant */
int skip_product_long(int last_most_significant, int width) {
   int test, temp_lms, offset;

   if (last_most_significant==0) {
      test=test_product(width, 0);
      if (test==0) return last_most_significant;
      if (test<0) {
         last_most_significant=overwrite_product(width, -test);
         if (last_most_significant==0) return last_most_significant;
      }
      else {
         last_most_significant=increment_product(width+1);
      }
   }

   offset=last_most_significant;
   while (offset>=0) {
      if (most_significant-width<offset) offset=most_significant-width;
      test=test_product(width+offset, offset);
      if (test==0) {
         offset--;
      }
      else if (test<0) {
         temp_lms=overwrite_product(width+offset, -test);
         if (temp_lms>last_most_significant) last_most_significant=temp_lms;
         if (temp_lms>offset) offset=temp_lms;
         else offset--;
      }
      else {
         temp_lms=increment_product(width+offset+1);
         if (temp_lms>last_most_significant) last_most_significant=temp_lms;
         offset=temp_lms;
      }
   }

   return last_most_significant;
}

int is_unique_product() {
   int i;
   
   for (i=0; i<num_unique; i++) {
      if (md(U1, unique_matrices[i]) > 4-epsilon) return 0;
   }

   return 1;
}

void calculate_product(int last_most_significant) {
   int i;

   if (last_most_significant==0) {
      U1=mm(U2,gate_list[product[0]]);
   }
   else if (last_most_significant==1) {
      U2=mm(U3,gate_list[product[1]]);
      U1=mm(U2,gate_list[product[0]]);
   }
   else if (last_most_significant==2) {
      U3=mm(U4,gate_list[product[2]]);
      U2=mm(U3,gate_list[product[1]]);
      U1=mm(U2,gate_list[product[0]]);
   }
   else if (last_most_significant==3) {
      U4=mm(U,gate_list[product[3]]);
      U3=mm(U4,gate_list[product[2]]);
      U2=mm(U3,gate_list[product[1]]);
      U1=mm(U2,gate_list[product[0]]);
   }
   else {
      mi(&U);
      for (i=most_significant; i>=4; i--) U=mm(U,gate_list[product[i]]);
      U4=mm(U,gate_list[product[3]]);
      U3=mm(U4,gate_list[product[2]]);
      U2=mm(U3,gate_list[product[1]]);
      U1=mm(U2,gate_list[product[0]]);
   }
}

void add_matrix_to_unique() {
   // TODO: get rid of unique matrix structure
   unique_matrices[num_unique]=U1;
   num_unique++;
}

void insert_product_in_list() {
   int i;

   unique_product_lists[free_list]=most_significant+1;
   for (i=0; i<=most_significant; i++) {
      unique_product_lists[free_list+1+i]=product[most_significant-i];
   }

   assert(unique_greater_than_product(free_list+1));

   free_list+=most_significant+2;
}

void insert_product_in_tree() {
   int i, current_position;
   
   current_position=0;

   current_position+=product[most_significant];
   for (i=1; i<=most_significant; i++) {
      current_position=product_check_tree[current_position];
      current_position+=product[most_significant-i];
   }

   product_check_tree[current_position]=free_node;
   product_check_tree[free_node]=current_position;
   free_node+=NG+1;

   current_position--;
   if (current_position > 0 && current_position%(NG+1)==0) {
      current_position--;
      for (i=1; i<=NG-1; i++)
         if (product_check_tree[current_position-NG+i]>0) {
            product_check_tree[current_position]=i;
            break;
         }
      current_position--;
   }
   while (current_position>0 && product_check_tree[current_position]==0) {
      product_check_tree[current_position]=-free_list;
      current_position--;
      if (current_position%(NG+1)==0) {
         current_position--;
         for (i=1; i<=NG-1; i++)
            if (product_check_tree[current_position-NG+i]>0) {
               product_check_tree[current_position]=i;
               break;
            }
         current_position--;
      }
   }
}

void add_matrix_to_structures() {
   insert_product_in_tree();
   insert_product_in_list();
   add_matrix_to_unique();
}
