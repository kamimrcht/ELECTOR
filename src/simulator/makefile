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


EXEC= simulator

all: $(EXEC)


simulator:   simulator.o
	$(CC) -o $@ $^ $(LDFLAGS)

simulator.o: simulator.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
