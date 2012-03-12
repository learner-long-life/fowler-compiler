#include "bin.h"

MatrixItem::MatrixItem(Matrix &_mtx, double _dist, int _sequence_index,
                       int _sequence_length, MatrixItem* _next) :
  mtx(_mtx), dist(_dist), sequence_index(_sequence_index),
  sequence_length(_sequence_length), next(_next) {}

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
  return bin;
}

void BinSet::insert(double dist, Matrix &m,
                    int sequence_index, int sequence_length,
                    bool side) {

  // Find place to insert item
  unsigned int bin = find(dist);
  MatrixItem** item = side ? &(left_bins[bin]) : &(right_bins[bin]);
  while (*item != NULL) { item = &((*item)->next); }

  // Add item to bin
  *item = new MatrixItem(m, dist, sequence_index,
                         sequence_length, NULL);
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
}

BinSet::~BinSet() {
  delete [] left_bins;
  delete [] right_bins;
}
