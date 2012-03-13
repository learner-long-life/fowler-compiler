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
  left_bins  = new MatrixItem*[bin_count];
  for (int i = 0; i < bin_count; i++) { left_bins[i] = NULL; }
  right_bins = new MatrixItem*[bin_count];
  for (int i = 0; i < bin_count; i++) { right_bins[i] = NULL; }
}

int BinSet::find(double dist) {
  int bin = (int)((dist - min) * bin_count / range);
  if (bin >= bin_count) { bin = bin_count - 1; }
  if (bin < 0) bin = 0;
  return bin;
}

void BinSet::insert(double dist, Matrix m,
                    int sequence_index, int sequence_length,
                    bool side) {

  // Find place to insert item
  unsigned int bin = find(dist);
  MatrixItem** item = side ? &(right_bins[bin]) : &(left_bins[bin]);
  while (*item != NULL) { item = &((*item)->next); }
#ifdef DEBUG_BINS
  cout << "Inserting item into bin " << bin << endl;
#endif

  // Add item to bin
  *item = new MatrixItem(m, dist, sequence_index,
                         sequence_length, NULL);
}

int BinSet::contains(Matrix m, double dist, double threshold, bool side,
                     double &acc) {
  // Determine number and location of bins to search in.
  int min_bin = find(dist - threshold);
  int max_bin = find(dist + threshold);
  int search_bin_count = max_bin - min_bin + 1;
#ifdef DEBUG_BINS
  cout << "Will search bins " << min_bin << " to " << max_bin << endl;
#endif

  // Initialize list of bin pointers.
  // They will keep track of our search position in each bin.
  MatrixItem **bins_to_search = new MatrixItem*[search_bin_count];
  int bin = 0;
  MatrixItem **haystack = side ? right_bins : left_bins;
  for (int i = min_bin; i <= max_bin; i++) {
    bins_to_search[bin++] = haystack[i];
  }

  // Search for the appropriate bin.
  int current_length = 1;
  bool items_left = true;
  int next_smallest_length;
  while (items_left) {
    items_left = false;
    next_smallest_length = numeric_limits<int>::max();
    for (int i = 0; i < search_bin_count; i++) {
      // See if we have an item to examine
      MatrixItem *item = bins_to_search[i];
      if (item == NULL) { continue; }
#ifdef DEBUG_BINS
      cout << "Examining matrix in bin " << (i+min_bin) << endl;
#endif
      items_left = true;

      // Skip over long sequences
      if (item->sequence_length < next_smallest_length) {
        next_smallest_length = item->sequence_length;
      }
      if (item->sequence_length > current_length) { continue; }

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
          delete [] bins_to_search;
          return item->sequence_index;
        }
      }

      // Since a result wasn't found, advance the pointer.
      bins_to_search[i] = item->next;
    }
    current_length = next_smallest_length;
  }
  delete [] bins_to_search;
  return -1;
}

void BinSet::destroy_bins(MatrixItem **bins) {
  for (int i = 0; i < bin_count; i++) {
    MatrixItem* item = bins[i];
    while (item != NULL) {
      MatrixItem* next = item->next;
      delete item;
      item = next;
    }
  }
  delete [] bins;
}

void BinSet::print_bin(MatrixItem **bin, FILE *out, void (*print_sequence)(FILE*,int)) {
  for (int i = 0; i < bin_count; i++) {
    fprintf(out, "===== Bin %d =====\n", i);
    MatrixItem* start = bin[i];
    while (start != NULL) {
      start->print(out, print_sequence);
      start = start->next;
    }
  }
}

void BinSet::print(FILE *out, void (*print_sequence)(FILE*,int)) {
  fprintf(out, "Left Bins\n");
  print_bin(left_bins, out, print_sequence);

  fprintf(out, "Right Bins\n");
  print_bin(right_bins, out, print_sequence);
}

BinSet::~BinSet() {
  destroy_bins(left_bins);
  destroy_bins(right_bins);
}
