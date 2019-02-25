all: ltq-nand

ltq-nand: ltq-nand.c
	gcc -o ltq-nand ltq-nand.c

clean:
	rm -f *.o ltq-nand
