#ifndef BIN_H
#define BIN_H

#include "matrix.h"

/**
 * A single matrix item.
 */
struct MatrixItem {
  Matrix mtx;         // Matrix to compare with
  double dist;        // Distance from the target matrix
  int sequence_index; // Index where the sequence that
                      // generated this matrix begins.
  int sequence_length; // Length of the sequence.
  struct MatrixItem* next;
};

/**
 * A set of bins, into which matrices are inserted by
 * distance.
 */
typedef struct {
  double min, range;
  int bin_count;
  struct MatrixItem **left_bins;
  struct MatrixItem **right_bins;
} BinSet;

/**
 * Creates a new BinSet.
 *  - min - the minimum distance value.
 *  - range - add to min to get the max distance value.
 *  - bin_count - the number of bins in this BinSet.
 */
BinSet* bin_create(double min, double range, int bin_count);

/**
 * Frees a BinSet.
 */
void bin_destroy(BinSet *bin_set);

#endif // BIN_H
