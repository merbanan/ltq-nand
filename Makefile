all: ltq-nand

ltq-nand: ltq-nand.c bch_lib.c
	gcc -O0 -g3 -o ltq-nand ltq-nand.c bch_lib.c

clean:
	rm -f *.o ltq-nand
