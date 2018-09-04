CC = g++
CFLAGS = -std=c++11 -O3 -Wall

binaries = main.o particle.o cycle.o

main: $(binaries)
	$(CC) $(CFLAGS) $(binaries) -o main -ljsoncpp

main.o: main.cc
	$(CC) $(CFLAGS) -c main.cc -o main.o

particle.o: particle.cc
	$(CC) $(CFLAGS) -c particle.cc -o particle.o

cycle.o: cycle.cc
	$(CC) $(CFLAGS) -c cycle.cc -o cycle.o

.PHONY: clean
clean:
	rm -rf main *.o