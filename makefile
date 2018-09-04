CC = g++
CFLAGS = -std=c++11 -O3 -Wall

binaries = main.o particle.o mesh.o

main: $(binaries)
	$(CC) $(CFLAGS) $(binaries) -o main -ljsoncpp

main.o: main.cc
	$(CC) $(CFLAGS) -c main.cc -o main.o

particle.o: particle.cc
	$(CC) $(CFLAGS) -c particle.cc -o particle.o

mesh.o: mesh.cc
	$(CC) $(CFLAGS) -c mesh.cc -o mesh.o

.PHONY: clean
clean:
	rm -rf main *.o