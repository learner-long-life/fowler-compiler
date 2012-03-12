#ifndef BIN_H
#define BIN_H

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
};

/**
 * A set of bins, into which matrices are inserted by
 * distance.
 */
class BinSet {
protected:
  double min, range;
  int bin_count;
  struct MatrixItem **left_bins;
  struct MatrixItem **right_bins;

  inline int find(double dist);

  /**
   * Destroys bin lists.
   */
  void destroy_bins(MatrixItem **bins);

public:

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
   */
  void insert(double dist, Matrix &m,
              int sequence_index, int sequence_length,
              bool side);
};

#endif // BIN_H
