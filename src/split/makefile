CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread -pipe -Wfatal-errors -fopenmp
LDFLAGS=-pthread -pipe -flto -fopenmp



EXEC= masterSplitter Donatello

all: $(EXEC)

sortPAF: sortPAF.o  utils.o mapping.o distances.o ssw_cpp.o ssw.o
	$(CC) -o $@ $^ $(LDFLAGS)

Donatello.o: Donatello.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

Donatello: Donatello.o
	$(CC) -o $@ $^ $(LDFLAGS)

masterSplitter.o: Master_Splitter.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

masterSplitter: masterSplitter.o
	$(CC) -o $@ $^ $(LDFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
