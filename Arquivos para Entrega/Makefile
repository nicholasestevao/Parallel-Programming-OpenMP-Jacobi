CC := gcc

ifeq ($(shell uname), Darwin)
	CC := gcc-13
endif

DBFLAGS := -Wall -g3

CFLAGS := -fopenmp -march=native -O3 $(DBFLAGS)

ifeq ($(OS),Windows_NT)
	OUT_EXT := .exe
else
	OUT_EXT :=
endif

all: seq par

seq: jacobiseq.c
	$(CC) $(CFLAGS) jacobiseq.c -o jacobiseq$(OUT_EXT)

par: jacobipar.c
	$(CC) $(CFLAGS) jacobipar.c -o jacobipar$(OUT_EXT)

runSeq: ./jacobiseq$(OUT_EXT)
	./jacobiseq$(OUT_EXT) $(ARGS)

runPar: ./jacobipar$(OUT_EXT)
	./jacobipar$(OUT_EXT) $(ARGS)

clean:
	rm -rf *.out *.exe *.txt

.PHONY: clean
