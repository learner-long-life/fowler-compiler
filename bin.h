#ifndef BIN_H
#define BIN_H

#include <stdio.h>
#include <map>
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

typedef std::map<int, MatrixItem**> BinMap;

/**
 * A set of bins, into which matrices are inserted by
 * distance.
 */
// TODO: remove the left bin because it's not used anymore.
// TODO: remove all but the longest and second-longest sequences to optimize
// memory use.
class BinSet {
protected:
  /**
   * Finds the index for the given distance.
   */
  inline int find(double dist);

  /**
   * Destroys the given bin.
   */
  void destroy_bin(MatrixItem** bin);

public:
  double min, range;
  int bin_count;
  BinMap bins;

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
   */
  void insert(double dist, Matrix m,
              int sequence_index, int sequence_length);

  /**
   * Determines if a matrix exists in the set that is within threshold of
   * the target matrix in the search, which is also within threshold of the
   * given matrix.
   *
   * If the matrix exists, the index of its gate sequence is returned.
   *
   * Otherwise, -1 is returned.
   */
  int contains(Matrix m, double dist, double threshold,
               double &acc);

  /**
   * Delete all sequences shorter than min_length, to speed up computation.
   */
  // TODO: compress the heap to improve cache locality.
  void delete_short_sequences(int min_length);
};

#endif // BIN_H
