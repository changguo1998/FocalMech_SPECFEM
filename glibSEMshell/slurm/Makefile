CFLAGS = -O3 -g -W

all : calcglib compress

compress : compress.c GlibIO.c Huffman.c
	gcc $(CFLAGS) compress.c -o ../semexec/bin/compress -lm -fopenmp

calcglib : calcglib.c
	gcc $(CFLAGS) calcglib.c -o ../semexec/bin/calcglib -fopenmp