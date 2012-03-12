#include "bin.h"

BinSet* bin_create(double min, double range, int bin_count) {
  BinSet* bin_set = new BinSet;
  if (bin_set == NULL) { return NULL; }

  bin_set->min        = min;
  bin_set->range      = range;
  bin_set->bin_count  = bin_count;
  bin_set->left_bins  = new MatrixItem*[bin_count];
  if (bin_set->left_bins == NULL) {
    delete bin_set;
    return NULL;
  }
  bin_set->right_bins  = new MatrixItem*[bin_count];
  if (bin_set->right_bins == NULL) {
    delete [] bin_set->left_bins;
    delete bin_set;
    return NULL;
  }
  return bin_set;
}

/**
 * Frees a BinSet.
 */
void bin_destroy(BinSet *bin_set) {
  delete [] bin_set->left_bins;
  delete [] bin_set->right_bins;
  delete bin_set;
}
