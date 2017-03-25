CC=g++
CFLAGS=-g -Wall -O3 -I/usr/local/inlcude

all: Mutations

Mutations: Main.o Pop.o
	$(CC) -L/usr/local/include -lgsl Main.o Pop.o -o Mutations
    
Main.o: Main.cpp Pop.h
	$(CC) $(CFLAGS) -c Main.cpp -o Main.o

Pop.o: Pop.cpp Pop.h
	$(CC) $(CFLAGS) -c Pop.cpp -o Pop.o
clean:
	rm *o Mutations