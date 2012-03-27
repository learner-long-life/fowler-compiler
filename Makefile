gate: complex.o matrix.o main07a.o
	cc -o gate complex.o matrix.o main07a.o

complex.o: complex.h complex.c
	cc -c complex.c

matrix.o: matrix.h matrix.c
	cc -c matrix.c

main07a.o: complex.h matrix.h main07a.c
	cc -c main07a.c
