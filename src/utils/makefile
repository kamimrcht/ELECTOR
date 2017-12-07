#~ CC=/usr/bin/g++
CC=g++
CFLAGS=  -Wall -Wextra  -Ofast -std=c++11 -march=native -pthread -pipe -Isparsepp
LDFLAGS=-pthread -Isparsepp


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=fa2fq


fa2fq:   fa2fq.o
	$(CC) -o $@ $^ $(LDFLAGS)

fa2fq.o: fatofq.cpp
	$(CC) -o $@ -c $< $(CFLAGS)
