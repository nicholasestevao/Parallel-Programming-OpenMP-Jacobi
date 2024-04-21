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

seq: jacobiseq.c
	$(CC) $(CFLAGS) jacobiseq.c -o jacobiseq$(OUT_EXT)

par: jacobipar.c
	$(CC) $(CFLAGS) jacobipar.c -o jacobipar$(OUT_EXT)

clean:
	rm -rf *.out *.exe

.PHONY: clean
