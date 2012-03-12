#!/bin/bash
# Runs unit tests on bin.
cd tests

# Test creating a structure.
valgrind ./bin_test 0 2>results.tmp

if [[ $? != 0 ]]; then
  echo "FAILED creating structure"
else
  if [[ -e create_structure_results.txt ]]; then
    sed "s/\s*==.*==\s*//" < results.tmp | diff - create_structure_results.txt
    if [[ $? == 1 ]]; then
      echo "FAILED creating structure"
    else
      echo "PASSED creating structure"
    fi
  else
    echo "Created structure results file:"
    sed "s/\s*==.*==\s*//" < results.tmp > create_structure_results.txt
    cat create_structure_results.txt
  fi
fi
