#ifndef BIN_H
#define BIN_H

#include <stdio.h>
#include "matrix.h"

/**
 * A single matrix item.
 */
class MatrixItem {
public:
  Matrix mtx;         // Matrix to compare with
  double dist;        // Distance from the target matrix
  int sequence_index; // Index where the sequence that
                      // generated this matrix begins.
  int sequence_length; // Length of the sequence.
  MatrixItem* next;

  MatrixItem(Matrix &_mtx, double _dist, int _sequence_index,
             int _sequence_length, MatrixItem* _next);

  void print(FILE *out, void (*print_sequence)(FILE*,int));
};

/**
 * A set of bins, into which matrices are inserted by
 * distance.
 */
class BinSet {
protected:
  inline int find(double dist);

  /**
   * Destroys bin lists.
   */
  void destroy_bins(MatrixItem **bins);

public:
  double min, range;
  int bin_count;
  struct MatrixItem **left_bins;
  struct MatrixItem **right_bins;

  /**
   * Creates a new BinSet.
   *  - min - the minimum distance value.
   *  - range - add to min to get the max distance value.
   *  - bin_count - the number of bins in this BinSet.
   */
  BinSet(double _min, double _range, int _bin_count);

  /**
   * Frees a BinSet.
   */
  ~BinSet();

  /**
   * Inserts a matrix into the bin set.
   * TODO: finish documentation here.
   *  - side - true if the matrix is on the right side,
   *           false on the left.
   */
  void insert(double dist, Matrix m,
              int sequence_index, int sequence_length,
              bool side);

  /**
   * Determines if a matrix exists in the set that is within threshold of
   * the target matrix in the search, which is also within threshold of the
   * given matrix.
   *
   * If the matrix exists, the index of its gate sequence is returned.
   *
   * Otherwise, -1 is returned.
   */
  int contains(Matrix m, double dist, double threshold, bool side,
               double &acc);

  void print_bin(MatrixItem **bin, FILE *out, void (*print_sequence)(FILE*, int));
  void print(FILE *out, void (*print_sequence)(FILE*, int));
};

#endif // BIN_H
