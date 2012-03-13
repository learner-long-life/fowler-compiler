#include "bin.h"
#include <limits>

#ifdef DEBUG_BINS
#include <iostream>
#endif

using namespace std;

MatrixItem::MatrixItem(Matrix &_mtx, double _dist, int _sequence_index,
                       int _sequence_length, MatrixItem* _next) :
  mtx(_mtx), dist(_dist), sequence_index(_sequence_index),
  sequence_length(_sequence_length), next(_next) {}

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

BinSet::BinSet(double _min, double _range, int _bin_count) {
  min        = _min;
  range      = _range;
  bin_count  = _bin_count;
}

int BinSet::find(double dist) {
  int bin = (int)((dist - min) * bin_count / range);
  if (bin >= bin_count) { bin = bin_count - 1; }
  if (bin < 0) bin = 0;
  return bin;
}

void BinSet::insert(double dist, Matrix m,
                    int sequence_index, int sequence_length) {
  // Get the bin list to insert the item into.
  MatrixItem** bin_list;
  if (bins.count(sequence_length) == 0) {
    bin_list = new MatrixItem*[bin_count];
    for (int i = 0; i < bin_count; i++) { bin_list[i] = 0; }
    bins[sequence_length] = bin_list;
  } else {
    bin_list = bins[sequence_length];
  }

  // Find place to insert item.
  unsigned int bin = find(dist);
  MatrixItem** item = &(bin_list[bin]);
  while (*item != NULL) { item = &((*item)->next); }
#ifdef DEBUG_BINS
  cout << "Inserting item into bin " << bin << endl;
#endif

  // Add item to bin
  *item = new MatrixItem(m, dist, sequence_index,
                         sequence_length, NULL);
}

void BinSet::delete_short_sequences(int min_length) {
  BinMap::iterator iter = bins.begin();
  while (iter != bins.end() && iter->first < min_length) {
    destroy_bin(iter->second);
    BinMap::iterator tmp = iter;
    iter++;
    bins.erase(tmp);
  }
}

int BinSet::contains(Matrix m, double dist, double threshold, double &acc) {
  // Determine number and location of bins to search in.
  int min_bin = find(dist - threshold);
  int max_bin = find(dist + threshold);
#ifdef DEBUG_BINS
  cout << "Will search bins " << min_bin << " to " << max_bin << endl;
#endif

  // Search over sequences by increasing length
  for (BinMap::iterator iter = bins.begin();
       iter != bins.end(); iter++) {
    // Search over all bins in the bin range
    for (int i = min_bin; i <= max_bin; i++) {
      MatrixItem *item = iter->second[i];
      while (item != NULL) {
        // Both matrices should have about the same distance
        // to the target matrix.  Otherwise, by the triangle
        // inequality, we know they won't be near each other.
        double dist_diff = dist - item->dist;
#ifdef DEBUG_BINS
          cout << "Matrix distance: " << dist_diff << endl;
#endif
        if (dist_diff >= -threshold && dist_diff <= threshold) {
          // Now, see if the matrices are within the threshold distance of
          // each other.  If they are, then one could connect their corresponding
          // sequences to obtain a result matrix.
          // Thus, the sequence for the discovered matrix should be returned.
          // TODO: optimize out the square roots
          double mtx_diff = md_tri(item->mtx, m);
#ifdef DEBUG_BINS
          cout << "Distance between matrices: " << mtx_diff << endl;
#endif
          if (mtx_diff <= threshold) {
#ifdef DEBUG_BINS
            cout << "MATRIX FOUND! Sequence " << item->sequence_index << endl;
#endif
            acc = mtx_diff;
            return item->sequence_index;
          }
        }
        item = item->next;
      }
    }
  }
  return -1;
}

void BinSet::destroy_bin(MatrixItem** bin) {
  for (int i = 0; i < bin_count; i++) {
    MatrixItem* item = bin[i];
    while (item != NULL) {
      MatrixItem* next = item->next;
      delete item;
      item = next;
    }
  }
  delete [] bin;
}

BinSet::~BinSet() {
  // Free all the bins in the bin map
  for (BinMap::iterator iter = bins.begin();
       iter != bins.end(); iter++) {
    destroy_bin(iter->second);
  }
}
