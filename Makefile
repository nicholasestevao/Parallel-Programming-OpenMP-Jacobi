CC := gcc

ifeq ($(shell uname), Darwin)
	CC := gcc-13
endif

DBFLAGS := -Wall -g3

CFLAGS := -fopenmp -march=native -O3 $(DBFLAGS)

ifeq ($(OS),Windows_NT)
	OUT_EXT := .exe
else
	OUT_EXT := .out
endif

all: seq par teste

seq: jacobiseq.c
	$(CC) $(CFLAGS) jacobiseq.c -o jacobiseq$(OUT_EXT)

par: jacobipar.c
	$(CC) $(CFLAGS) jacobipar.c -o jacobipar$(OUT_EXT)

teste: seq par teste.c
	$(CC) $(CFLAGS) teste.c -o teste$(OUT_EXT)

run: ./teste$(OUT_EXT) ./jacobiseq$(OUT_EXT) ./jacobipar$(OUT_EXT)
	./teste$(OUT_EXT)

clean:
	rm -rf *.out *.exe *.txt

.PHONY: clean
