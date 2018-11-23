CC = g++
CFLAGS = -std=c++11 -O3 -Wall

binaries = main.o functions.o FFT.o


main: $(binaries)
	$(CC) $(CFLAGS) $(binaries) -o main -ljsoncpp

main.o: main.cc functions.h
	$(CC) $(CFLAGS) -c main.cc -o main.o

functions.o: functions.cc functions.h
	$(CC) $(CFLAGS) -c functions.cc -o functions.o

FFT.o: FFT.cc functions.h
	$(CC) $(CFLAGS) -c FFT.cc -o FFT.o


.PHONY: clean
clean:
	rm -rf main *.o