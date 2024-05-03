CC := gcc

ifeq ($(shell uname), Darwin)
	CC := gcc-13
endif

DBFLAGS := -Wall -g3

CFLAGS := -fopenmp -march=native -O3 $(DBFLAGS)

ifeq ($(OS),Windows_NT)
	OUT_EXT :=
else
	OUT_EXT := .out
endif

all: seq par

seq: jacobi_sequencial.c
	$(CC) $(CFLAGS) jacobi_sequencial.c -o jacobi_sequencial$(OUT_EXT)

par: jacobi_paralelo.c
	$(CC) $(CFLAGS) jacobi_paralelo.c -o jacobi_paralelo$(OUT_EXT)

clean:
	rm -rf *.out *.exe

.PHONY: clean
