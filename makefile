CC = g++
CFLAGS = -std=c++11 -O3 -Wall

binaries = main.o functions.o


main: $(binaries)
	$(CC) $(CFLAGS) $(binaries) -o main -ljsoncpp

main.o: main.cc init.h
	$(CC) $(CFLAGS) -c main.cc -o main.o

functions.o: functions.cc init.h
	$(CC) $(CFLAGS) -c functions.cc -o functions.o


.PHONY: clean
clean:
	rm -rf main *.o