CC=mpicc

CFLAGS=-c -fPIC -Wall



all: a clean



a: main.o matrix_print.o matrix_init.o matrix_inverse.o norm.o 

	$(CC) main.o matrix_init.o matrix_print.o matrix_inverse.o norm.o  -o a -lm



main.o: main.c

	$(CC) $(CFLAGS) main.c



matrix_print.o: matrix_print.c

	$(CC) $(CFLAGS) matrix_print.c



matrix_init.o: matrix_init.c

	$(CC) $(CFLAGS) matrix_init.c

matrix_norm.o: norm.c
	$(CC) $(CFLAGS) norm.c

matrix_inverse.o: matrix_inverse.c
	$(CC) $(CFLAGS) matrix_print.h matrix_print.c matrix_inverse.c
clean:
	rm -rf *.o
