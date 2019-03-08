# ltq-nand
Reed Solomon encoder for the GRX350/GRX550 nand controller

Only BCH 512 7 bytes mode has been verified.

ltq-nand version 1.0
Usage: ltq-nand [OPTION]...
	 -i input file
	 -o output file
	 -b block size
	 -p page size
	 -s spare area/oob area size
	 -e ecc mode
		 0 0x00 fill
		 1 0xFF fill
		 2 Reed Solomon 3 byte
		 3 Reed Solomon 4 byte
		 4 BCH 512 7 bytes
		 5 BCH 512 13 bytes
Default command is:
./ltq-nand -i infile -o outfile -b 131072 -p 2048 -s 64 -m 4
