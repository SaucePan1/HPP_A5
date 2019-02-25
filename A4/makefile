CFLAGS= -Wall -O3 -ffast-math -march=native -msse4.2
LDFLAGS=-lm

galsim:sr_tree_bernat10.c
	gcc $(CFLAGS) $(INCLUDES) -c sr_tree_bernat10.c
	gcc -o galsim sr_tree_bernat10.o $(LDFLAGS)

make clean:
	rm -f ./galsim *.o
