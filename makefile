CC = g++
CFLAGS = -std=c++11 -O3

binaries = main.o particle.o functions.o


main: $(binaries)
	$(CC) $(CFLAGS) $(binaries) -o main

main.o: cpp-version/main.cc
	$(CC) $(CFLAGS) -c cpp-version/main.cc -o main.o

particle.o: cpp-version/particle.cc
	$(CC) $(CFLAGS) -c cpp-version/particle.cc -o particle.o

functions.o: cpp-version/functions.cc
	$(CC) $(CFLAGS) -c cpp-version/functions.cc -o functions.o


.PHONY: clean
clean:
	rm -rf main *.o