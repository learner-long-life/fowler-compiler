#include "bin.h"
#include <limits>

#ifdef DEBUG_BINS
#include <iostream>
#endif

using namespace std;

#ifdef COUNT_SEARCHES
unsigned long search_counter = 0, positive_counter = 0;
#endif

MatrixItem::MatrixItem(Matrix &_mtx, double _dist, int _sequence_index,
                       int _sequence_length) :
  mtx(_mtx), dist(_dist), sequence_index(_sequence_index),
  sequence_length(_sequence_length) {}

void MatrixItem::print(FILE *out, void (*print_sequence)(FILE*,int)) {
  fprintf(out, "-----------------\n");

  // Print sequence
  fprintf(out, "Sequence: ");
  print_sequence(out, sequence_index);

  // Print distance
  fprintf(out, "Distance: %.10f\n", dist);

  // Print matrix
  fprintf(out, "Matrix:\n");
  pm(out, mtx);
}

BinSet::BinSet() {}

void BinSet::insert(double dist, Matrix m,
                    int sequence_index, int sequence_length) {
  DistMap &mtx_set = bins[sequence_length];
  // TODO: does MatrixItem need a dist property anymore?
  mtx_set.insert(DistMapPair(dist,
                 MatrixItem(m, dist, sequence_index, sequence_length)));
}

void BinSet::delete_short_sequences(int min_length) {
  BinMap::iterator iter = bins.begin();
  while (iter != bins.end() && iter->first < min_length) {
    BinMap::iterator tmp = iter;
    iter++;
    bins.erase(tmp);
  }
}

int BinSet::contains(Matrix m, double dist, double threshold, double &acc) {
#ifdef COUNT_SEARCHES
    int ret_val = -1;
#endif

  // Search over sequences by increasing length
  for (BinMap::iterator iter = bins.begin();
       iter != bins.end(); iter++) {
    DistMap &mtx_set = iter->second;
    // Find lower bound
    DistMapIter iter = mtx_set.lower_bound(dist - threshold);

    // Find upper bound
    // TODO: what does this return if there is no element?  I assume end()
    DistMapIter end = mtx_set.upper_bound(dist + threshold);

    // Search over all bins in the bin range
    for (; iter != end; iter++) {
      MatrixItem &item = iter->second;

      // Both matrices should have about the same distance
      // to the target matrix.  Otherwise, by the triangle
      // inequality, we know they won't be near each other.
#ifdef DEBUG_BINS
      cout << "Matrix distance: " << dist_diff << endl;
#endif
      // Now, see if the matrices are within the threshold distance of
      // each other.  If they are, then one could connect their corresponding
      // sequences to obtain a result matrix.
      // Thus, the sequence for the discovered matrix should be returned.
      // TODO: optimize out the square roots
      double mtx_diff = md_tri(item.mtx, m);

#ifdef DEBUG_BINS
        cout << "Distance between matrices: " << mtx_diff << endl;
#endif

#ifdef COUNT_SEARCHES
      ret_val = item.sequence_index;
      search_counter++;
#endif

      if (mtx_diff <= threshold) {
#ifdef DEBUG_BINS
        cout << "MATRIX FOUND! Sequence " << item.sequence_index << endl;
#endif
        acc = mtx_diff;
#ifdef COUNT_SEARCHES
        ret_val = item.sequence_index;
        positive_counter++;
#else
        return item.sequence_index;
#endif
      }
    }
  }

#ifdef COUNT_SEARCHES
  return ret_val;
#else
  return -1;
#endif
}

