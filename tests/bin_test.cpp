#include <cstdlib>
#include <cmath>
#include <iostream>
#include "../bin.h"

using namespace std;

/**
 * Test creating and destroying a BinSet.
 */
int test_creating_bin_set() {
  BinSet *set = new BinSet(0, 4, 10000);
  delete set;
  return 0;
}

/**
 * Test inserting an item into a BinSet, and retrieving it.
 */
int test_inserting() {
  BinSet set(0, 4, 10000);
  Matrix m;
  double sqrt2o2 = sqrt(2)/2;
  m.z11.x = sqrt2o2;
  m.z11.y = 0;
  m.z12.x = sqrt2o2;
  m.z12.y = 0;
  m.z21.x = sqrt2o2;
  m.z21.y = 0;
  m.z22.x = -sqrt2o2;
  m.z22.y = 0;
  set.insert(2, m, 10, 20, false);
  return set.contains(m, 2, .01, false) == 10 ? 0 : 1;
}

int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 2) {
    cerr << "Missing test number." << endl;
    return -1;
  }

  // Run test
  int test_number = atoi(argv[1]);
  switch (test_number) {
    case  0: return test_creating_bin_set();
    case  1: return test_inserting();
    default: cerr << "Invalid test number." << endl; return -2;
  }

  return 0;
}
