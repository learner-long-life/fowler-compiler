#!/bin/bash
# Runs unit tests on bin.
cd tests

# Runs a test.  Provide this function with the test number.
run_valgrind_test() {
  # Test creating a structure.
  valgrind ./bin_test "$1" 2>results.tmp

  if [[ $? != 0 ]]; then
    echo "FAILED $2"
  else
    if [[ -e "$1_stderr.txt" ]]; then
      sed "s/\s*==.*==\s*//" < results.tmp | diff - "$1_stderr.txt"
      if [[ $? == 1 ]]; then
        echo "FAILED $2"
      else
        echo "PASSED $2"
      fi
    else
      echo "Created structure results file:"
      sed "s/\s*==.*==\s*//" < results.tmp > "$1_stderr.txt"
      cat "$1_stderr.txt"
    fi
  fi
}

run_valgrind_test 0 "creating structure"
run_valgrind_test 1 "inserting and finding in structure"
