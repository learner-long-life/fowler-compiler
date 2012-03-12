#include <iostream>
#include <cstdlib>
#include "../bin.h"

using namespace std;

/**
 * Test creating and destroying a BinSet.
 */
void test_creating_bin_set() {
  BinSet *set = new BinSet(0, 4, 10000);
  delete set;
}

int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 2) {
    cerr << "Missing test number." << endl;
    return 1;
  }

  // Run test
  int test_number = atoi(argv[1]);
  switch (test_number) {
    case 0: test_creating_bin_set(); break;
    default: cerr << "Invalid test number." << endl; return 2;
  }

  return 0;
}
