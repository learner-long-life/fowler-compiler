gate: complex.o matrix.o main07a.o
	g++ -lm -o gate complex.o matrix.o main07a.o

complex.o: complex.h complex.c
	g++ -g3 -c complex.c

matrix.o: matrix.h matrix.c
	g++ -g3 -c matrix.c

main07a.o: complex.h matrix.h main07a.c
	g++ -g3 -c main07a.c

clean:
	rm -f gate complex.o matrix.o main07a.o *.tmp

gate_generator.o: gate_generator.c
	g++ -c gate_generator.c

gate_generator: gate_generator.o complex.o matrix.o
	g++ -lm -o gate_generator complex.o matrix.o gate_generator.o

gates_37bit.mif: gate_generator
	./gate_generator

matrix_benchmark.o: matrix_benchmark.c
	g++ -c matrix_benchmark.c

matrix_benchmark: matrix_benchmark.o complex.o matrix.o
	g++ -lm -lrt -o matrix_benchmark complex.o matrix.o matrix_benchmark.o

# Benchmarking
gate_bench: complex.o matrix.o main07a_bench.o
	g++ -lm -lrt -o gate_bench complex.o matrix.o main07a_bench.o

main07a_bench.o: complex.h matrix.h main07a.c
	g++ -g3 -DBENCHMARK -c main07a.c -o main07a_bench.o

# Matrix distance calculating
gate_dist: complex.o matrix.o main07a_dist.o
	g++ -lm -lrt -o gate_dist complex.o matrix.o main07a_dist.o

main07a_dist.o: complex.h matrix.h main07a.c
	g++ -g3 -DDISTANCES -c main07a.c -o main07a_dist.o

# Meet-in-the-middle binning structure
bin.o: bin.cpp bin.h
	g++ -g3 -c bin.cpp -o bin.o

# Test binning structure
tests/bin_debug.o: bin.cpp bin.h
	g++ -g3 -DDEBUG_BINS -c bin.cpp -o bin.o

tests/bin_test.o: tests/bin_test.cpp
	g++ -g3 -c tests/bin_test.cpp -o tests/bin_test.o

tests/bin_test: tests/bin_test.o tests/bin_debug.o matrix.o complex.o
	g++ -lm -lrt -o tests/bin_test tests/bin_test.o tests/bin_debug.o matrix.o complex.o

test: tests/bin_test
	tests/tests.sh

# Use binning structure with profiling
gate_bin_fast: complex.o matrix.o main07a_bin_fast.o bin.o
	g++ -lm -lrt -o gate_bin_fast complex.o matrix.o main07a_bin_fast.o bin.o

main07a_bin_fast.o: complex.h matrix.h main07a.c
	g++ -pg -g3 -DBENCHMARK -DBIN -c main07a.c -o main07a_bin_fast.o

# Use binning structure with profiling
gate_bin: complex_profile.o matrix_profile.o main07a_bin.o bin.o
	g++ -pg -lm -lrt -o gate_bin complex_profile.o matrix_profile.o main07a_bin.o bin.o

main07a_bin.o: complex.h matrix.h main07a.c
	g++ -pg -g3 -DBENCHMARK -DBIN -DFIRST_STAGE_ONLY -c main07a.c -o main07a_bin.o

# Profiling
gate_profile: complex_profile.o matrix_profile.o main07a_profile.o
	g++ -pg -lm -o gate_profile complex_profile.o matrix_profile.o main07a_profile.o

gate_profile_long: complex_profile.o matrix_profile.o main07a.o
	g++ -pg -lm -o gate_profile_long complex_profile.o matrix_profile.o main07a.o

complex_profile.o: complex.h complex.c
	g++ -pg -g3 -c complex.c -o complex_profile.o

matrix_profile.o: matrix.h matrix.c
	g++ -pg -g3 -c matrix.c -o matrix_profile.o

main07a_profile.o: complex.h matrix.h main07a.c
	g++ -pg -g3 -DFIRST_STAGE_ONLY -c main07a.c -o main07a_profile.o

# Use SIMD
gate_simd: complex_simd.o matrix.o main07a_simd.o
	g++ -lm -lrt -o gate_simd complex_simd.o matrix.o main07a_simd.o

main07a_simd.o: complex.h matrix.h main07a.c
	g++ -g3 -DBENCHMARK -c main07a.c -o main07a_simd.o

complex_simd.o: complex.h complex.c
	g++ -g3 -DCOMPLEX_SIMD -c complex.c -o complex_simd.o

# Use inlining
gate_inline: main07a_inline.o
	g++ -lm -lrt -o gate_inline main07a_inline.o

main07a_inline.o: main07a.c
	g++ -g3 -DUSE_INLINE -DBENCHMARK -c main07a.c -o main07a_inline.o

# Use product lookup tree
gate_plt: main07a_plt.o complex.o matrix.o
	g++ -lm -lrt -o gate_plt main07a_plt.o complex.o matrix.o

main07a_plt.o: main07a.c
	g++ -g3 -DPRODUCT_LOOKUP_TREE -DBENCHMARK -c main07a.c -o main07a_plt.o

# Use binning with alternate distance reference point.
gate_bin_2: main07a_bin_2.o complex.o matrix.o bin.o
	g++ -lm -lrt -o gate_bin_2 main07a_bin_2.o complex.o matrix.o bin.o

main07a_bin_2.o: main07a.c
	g++ -g3 -DOTHER_DIST -DBIN -DBENCHMARK -c main07a.c -o main07a_bin_2.o

# Use Aram Harrow's modified Pauli basis
gate_pb: main07a_pb.o complex.o matrix_pb.o
	g++ -lm -lrt -o gate_pb main07a_pb.o matrix_pb.o complex.o

main07a_pb.o: main07a.c
	g++ -g3 -DPAULI_BASIS -c main07a.c -o main07a_pb.o

matrix_pb.o: matrix.c
	g++ -g3 -DPAULI_BASIS -c matrix.c -o matrix_pb.o
