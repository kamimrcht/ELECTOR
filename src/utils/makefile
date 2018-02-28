CC=g++
CFLAGS=  -Wall -Wextra  -Ofast -std=c++11 -march=native -pthread -pipe -Isparsepp
LDFLAGS=-pthread -Isparsepp




EXEC=fq2fa


fq2fa:   fq2fa.o
	$(CC) -o $@ $^ $(LDFLAGS)

fq2fa.o: fq2fa.cpp
	$(CC) -o $@ -c $< $(CFLAGS)
