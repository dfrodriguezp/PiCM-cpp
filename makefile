CC = g++
CFLAGS = -std=c++11 -O3 -Wall

binaries = main.o particle.o functions.o FFT.o


main: $(binaries)
	$(CC) $(CFLAGS) $(binaries) -o main -lhdf5_serial

main.o: main.cc
	$(CC) $(CFLAGS) -c main.cc -o main.o

particle.o: particle.cc
	$(CC) $(CFLAGS) -c particle.cc -o particle.o

functions.o: functions.cc
	$(CC) $(CFLAGS) -c functions.cc -o functions.o

FFT.o: FFT.cc
	$(CC) $(CFLAGS) -c FFT.cc -o FFT.o

.PHONY: clean
clean:
	rm -rf main *.o