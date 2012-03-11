gate: complex.o matrix.o main07a.o
	g++ -lm -o gate complex.o matrix.o main07a.o

complex.o: complex.h complex.c
	g++ -g3 -c complex.c

matrix.o: matrix.h matrix.c
	g++ -g3 -c matrix.c

main07a.o: complex.h matrix.h main07a.c
	g++ -g3 -c main07a.c

clean:
	rm -f gate complex.o matrix.o main07a.o

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

# Matrix distance binning
#gate_dist: complex.o matrix.o main07a_dist.o
#	g++ -lm -lrt -o gate_dist complex.o matrix.o main07a_dist.o
#
#main07a_dist.o: complex.h matrix.h main07a.c
#	g++ -g3 -DDISTANCES -c main07a.c -o main07a_dist.o