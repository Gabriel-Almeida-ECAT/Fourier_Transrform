CC = gcc
CFLAGS = -Wall -O2 -w
LDFLAGS = -lm

all: fourier_algorithms

fourier_algorithms: fourier_algorithms.c
	$(CC) $(CFLAGS) -o fourier_algorithms fourier_algorithms.c $(LDFLAGS)

clean:
	rm -f fourier_algorithms
