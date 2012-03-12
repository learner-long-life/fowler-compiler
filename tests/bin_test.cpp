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

/**
 * Creates a pi*n gate.
 */
Matrix get_pi_gate(double n) {
  Matrix G;
  double x;
  G.z11.x=1;
  x = acos(-1)*n;
  G.z22.x=cos(x);
  G.z22.y=sin(x);
  return G;
}

/**
 * Test inserting several items into the same bin.
 */
int test_inserting_same_bin() {
  BinSet set(0, 4, 10000);
  Matrix pi_over_6_gate = get_pi_gate(1./6);
  Matrix pi_over_2_gate = get_pi_gate(1./2);
  set.insert(2, pi_over_2_gate, 10, 20, false);
  set.insert(2, pi_over_6_gate, 10, 20, false);
  return set.contains(pi_over_6_gate, 2, .01, false) == 10 ? 0 : 1;
}

/**
 * Test inserting into the right bin.
 */
int test_inserting_right_bin() {
  BinSet set(0, 4, 10000);
  Matrix pi_over_2_gate = get_pi_gate(1./2);
  set.insert(2, pi_over_2_gate, 10, 20, true);
  return set.contains(pi_over_2_gate, 2, .01, true) == 10 ? 0 : 1;
}

/**
 * Ensure that the bin with no elements is empty.
 */
int test_inserting_opposite() {
  BinSet set(0, 4, 10000);
  Matrix pi_over_2_gate = get_pi_gate(1./2);
  set.insert(2, pi_over_2_gate, 10, 20, true);
  return set.contains(pi_over_2_gate, 2, .01, false) == -1 ? 0 : 1;
}

/**
 * Test searching for the wrong matrix.
 */
int test_searching_wrong_matrix() {
  BinSet set(0, 4, 10000);
  Matrix pi_over_6_gate = get_pi_gate(1./6);
  Matrix pi_over_2_gate = get_pi_gate(1./2);
  set.insert(2, pi_over_2_gate, 10, 20, false);
  return set.contains(pi_over_6_gate, 2, .01, false) == -1 ? 0 : 1;
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
    case  2: return test_inserting_same_bin();
    case  3: return test_inserting_right_bin();
    case  4: return test_inserting_opposite();
    case  5: return test_searching_wrong_matrix();
    default: cerr << "Invalid test number." << endl; return -2;
  }

  return 0;
}
