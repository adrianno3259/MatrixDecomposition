CC=gcc
LIB=-lm
all: main.c matrix.c matrix.h
	$(CC) *.c -o main $(LIB)
