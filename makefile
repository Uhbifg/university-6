CC=mpicc

CFLAGS=-c -fPIC -Wall -DDEBUG



all: a clean



a: main.o matrix_print.o matrix_init.o 

	$(CC) main.o matrix_init.o matrix_print.o -o a -lm



main.o: main.c

	$(CC) $(CFLAGS) main.c



matrix_print.o: matrix_print.c

	$(CC) $(CFLAGS) matrix_print.c



matrix_init.o: matrix_init.c

	$(CC) $(CFLAGS) matrix_init.c




clean:
	rm -rf *.o
