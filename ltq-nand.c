/*
 * Reed Solomon encoder for the GRX350/GRX550 nand controller
 *
 * Copyright (C) 2018 iopsys Software Solutions AB. All rights reserved.
 *
 * Author: Benjamin Larsson <benjamin.larsson@iopsys.eu>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * version 2 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_BLOCK_SIZE 524288
#define MAX_OOB_BLOCK_SIZE 540672

#define ONES 1
#define ZEROS 0
#define ECC4 2
#define ECC3 3

const char *oob_mode_string[4] = { "Zeros", "Ones", "ECC4", "ECC3" };

unsigned char block_buffer[MAX_BLOCK_SIZE];
unsigned char out_buffer[MAX_OOB_BLOCK_SIZE];


unsigned char g_num2alpha[256] = {0,0,1,99,2,198,100,106,3,205,199,188,101,126,107,42,4,141,206,
    78,200,212,189,225,102,221,127,49,108,32,43,243,5,87,142,232,207,
    172,79,131,201,217,213,65,190,148,226,180,103,39,222,240,128,177,
    50,53,109,69,33,18,44,13,244,56,6,155,88,26,143,121,233,112,208,
    194,173,168,80,117,132,72,202,252,218,138,214,84,66,36,191,152,149,
    249,227,94,181,21,104,97,40,186,223,76,241,47,129,230,178,63,51,
    238,54,16,110,24,70,166,34,136,19,247,45,184,14,61,245,164,57,59,7,
    158,156,157,89,159,27,8,144,9,122,28,234,160,113,90,209,29,195,123,
    174,10,169,145,81,91,118,114,133,161,73,235,203,124,253,196,219,30,
    139,210,215,146,85,170,67,11,37,175,192,115,153,119,150,92,250,82,
    228,236,95,74,182,162,22,134,105,197,98,254,41,125,187,204,224,211,
    77,140,242,31,48,220,130,171,231,86,179,147,64,216,52,176,239,38,
    55,12,17,68,111,120,25,154,71,116,167,193,35,83,137,251,20,93,248,
    151,46,75,185,96,15,237,62,229,246,135,165,23,58,163,60,183};

unsigned char g_alpha2num[256] = {1,2,4,8,16,32,64,128,135,137,149,173,221,61,122,244,111,222,59,118,
    236,95,190,251,113,226,67,134,139,145,165,205,29,58,116,232,87,174,
    219,49,98,196,15,30,60,120,240,103,206,27,54,108,216,55,110,220,63,
    126,252,127,254,123,246,107,214,43,86,172,223,57,114,228,79,158,187,
    241,101,202,19,38,76,152,183,233,85,170,211,33,66,132,143,153,181,
    237,93,186,243,97,194,3,6,12,24,48,96,192,7,14,28,56,112,224,71,142,
    155,177,229,77,154,179,225,69,138,147,161,197,13,26,52,104,208,39,
    78,156,191,249,117,234,83,166,203,17,34,68,136,151,169,213,45,90,
    180,239,89,178,227,65,130,131,129,133,141,157,189,253,125,250,115,
    230,75,150,171,209,37,74,148,175,217,53,106,212,47,94,188,255,121,
    242,99,198,11,22,44,88,176,231,73,146,163,193,5,10,20,40,80,160,
    199,9,18,36,72,144,167,201,21,42,84,168,215,41,82,164,207,25,50,
    100,200,23,46,92,184,247,105,210,35,70,140,159,185,245,109,218,51,
    102,204,31,62,124,248,119,238,91,182,235,81,162,195,1};

/* Galois Field mul */
unsigned char g_mul(unsigned char arg1, unsigned char arg2) {
    unsigned char s=0;
    if ((arg1== 0) || (arg2==0)) {
        return 0;
    };
    s = (s + g_num2alpha[arg1]) % 255;
    s = (s + g_num2alpha[arg2]) % 255;
    return g_alpha2num[s];
}

/* Galois Field add */
unsigned char g_add(unsigned char  arg3, unsigned char  arg4) {
    unsigned char  s=0;
    s = s ^ arg3;
    s = s ^ arg4;
    return s;
};

void reed_solomon_128bytes_ecc(unsigned char *data_bytes_partial, unsigned char *s, int strength) {
    unsigned char g[4] = {0};
    unsigned char temp[4] = {0};
    unsigned char degree;
    unsigned char bytes;
    unsigned char y;
    int i;
    unsigned char t[3] = {0};
    
    if (strength == ECC4) {
        g[3] = 205;
        g[2] =  63;
        g[1] =  92;
        g[0] =  32;
        degree = 4;
        bytes = 128;

        for (i=0;i<=bytes-1;i++) {
            y     = g_add( s[3], data_bytes_partial[i] );
            temp[1]  = g_mul( y, g[3] );
            temp[2]  = g_mul( y, g[2]);
            temp[3]  = g_mul( y, g[1] );
            s[3]  = g_add( s[2], temp[1] );
            s[2]  = g_add( s[1], temp[2] );
            s[1]  = g_add( s[0], temp[3] );
            s[0]  = g_mul( y, g[0] );
        };
    }
    
    /* 3 bytes mode */
    if (strength == ECC3) {

        g[3] = 0;
        g[2] = g_add(g_add(2,4), 8);
        g[1] = g_add(g_mul(2,4), g_add(g_mul(4,8), g_mul(2,8)));
        g[0] = g_mul(g_mul(2,4),8);
        bytes = 128;
        for (i=0;i<=bytes-1;i++) {
            y   = g_add( s[2], data_bytes_partial[i] );
            temp[2]  = g_mul( y, g[2]);
            temp[3]  = g_mul( y, g[1] );
            s[2]  = g_add( s[1], temp[2] );
            s[1]  = g_add( s[0], temp[3] );
            s[0]  = g_mul( y, g[0] );
        };
    }
    
    return;
};

/*
for(i=0;i<mtd->writesize/128;i++){
    s = reed_solomn_128bytes_ecc(buf+i*128);
    chip->oob_poi[i*4]=s[3];
    chip->oob_poi[i*4+1]=s[2];
    chip->oob_poi[i*4+2]=s[1];
    chip->oob_poi[i*4+3]=s[0];
   
    
      for(i=0;i<mtd->writesize/128;i++){
             s = reed_solomn_128bytes_ecc(buf+i*128);
             chip->oob_poi[i*3]=s[2];
             chip->oob_poi[i*3+1]=s[1];
             chip->oob_poi[i*3+2]=s[0];
         }

    
};
*/


void print_usage() {
    printf("ltq-nand version 1.0\n");     
    printf("Usage: ltq-nand [OPTION]...\n");
    printf("\t -i input file\n");
    printf("\t -o output file\n");
    printf("\t -b block size\n");
    printf("\t -p page size\n");
    printf("\t -1 fill spare area/OOB with 0xFF\n");
    printf("\t -0 fill spare area/OOB with 0x00\n");
}

int main(int argc,char** argv)
{
    int option = 0;

    int block_size = 131072;
    int blocks;
    int b_idx, p_idx, ei_idx;
    int page_size = 2048;
    int pages_in_block = 0;
    int padd_zero = 0;
    int padd_ones = 0;
    int ecc_strength = 4;
    int oob_mode = ECC4;
    int oob_area_size = 64;
    const char* input_file = NULL;
    const char* output_file = NULL;
    FILE *in_file, *out_file;
    struct stat st;

    while ((option = getopt(argc, argv,"e:i:o:b:p:10")) != -1) {
        switch (option) {
            case 'i' : input_file = optarg; 
                break;
            case 'o' : output_file = optarg;
                break;
            case 'b' : block_size = atoi(optarg);
                break;
            case '1' : padd_ones = 1;
                break;
            case '0' : padd_zero = 1;
                break;
            case 'p' : page_size = atoi(optarg);
                break;
            case 'e' : ecc_strength = atoi(optarg);
                break;
            default: print_usage(); 
                 exit(EXIT_FAILURE);
        }
    }
    if (argc < 4) {
        print_usage();
        goto exit;
    }
    
    /* Validate input */
    if (block_size == 0) {
        printf("Block size can not be 0\n");
        goto exit;
    }
    if (block_size > MAX_BLOCK_SIZE) {
        printf("Block size to large\n");
    }
    
    if (page_size == 0) {
        printf("Page size can not be 0\n");
        goto exit;
    }
    if (page_size % 32) {
        printf("Page size not even multiple of 32\n");
        goto exit;
    } 
    
    
    if ((ecc_strength != 3) && (ecc_strength != 4)) {
        printf("-e only ecc strength of 3 or 4 bytes supported.\n");
        goto exit;
    }
    
    if (input_file == NULL) {
        printf("-i input file name is missing\n");
        goto exit;
    }
    if (output_file == NULL) {
        printf("-o output file name is missing\n");
        goto exit;
    }
    
    in_file = fopen(input_file,"r");
    if (in_file == NULL) {
        printf("In file open error\n");
        goto exit;
    }
    out_file = fopen(output_file,"wb");
    if (out_file == NULL) {
        printf("Out file open error\n");
        goto exit;
    }
    
    /* Validate that the input file size is a multiple of the block size */
    if (stat(input_file, &st) != 0) {
        printf("File size error\n");
        goto exit;
    }
    
    if ((st.st_size % block_size) != 0) {
        printf("File size not multiple of block size\n");
        goto exit;
    }
    
    /* configure oob area data */
    if (padd_zero)
        oob_mode = ZEROS;
    if (padd_ones)
        oob_mode = ONES;
    if (ecc_strength == 3)
        oob_mode = ECC3;
    
    
    blocks = st.st_size / block_size;
    oob_area_size = page_size / 32;
    pages_in_block = block_size / page_size;
    
    printf("Input file %s\n", input_file);
    printf("Size of input file: %d\n", (int)st.st_size);
    printf("Block size: %d\n", block_size);
    printf("Input file blocks %d\n", blocks);
    printf("Page size: %d\n", page_size);
    printf("Pages in block: %d\n", pages_in_block);
    printf("ECC strength: %d\n", ecc_strength);
    printf("OOB area generation: %s\n", oob_mode_string[oob_mode]);
    printf("OOB area size %d\n", oob_area_size);
    
    for (b_idx = 0 ; b_idx < blocks ; b_idx++ ) {
        //printf("Block %d\n", b_idx);
        fread(block_buffer, block_size, 1, in_file);
        for (p_idx = 0 ; p_idx < pages_in_block ; p_idx++) {
            //printf("\tPage %d\n", p_idx);
            /* copy the initial page */
            memcpy(&out_buffer[p_idx*(page_size+oob_area_size)], &block_buffer[p_idx*page_size], page_size);
            /* generate oob area data */
            //printf("offset %d 0x%x\n", (p_idx+1)*page_size, (p_idx+1)*page_size);
            switch (oob_mode) {
                case ZEROS: memset(&out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size], 0, oob_area_size); break;
                case ONES:  memset(&out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size], 0xFF, oob_area_size); break;
                case ECC4: {
                    memset(&out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size], 0xFF, oob_area_size);
                    for (ei_idx = 0 ; ei_idx < page_size/128 ; ei_idx++) {
                        unsigned char s[4] = {0};
                        reed_solomon_128bytes_ecc(&block_buffer[ p_idx*page_size + ei_idx*128], s, ECC4);
                        out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size + ei_idx*4  ] = s[3];
                        out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size + ei_idx*4+1] = s[2];
                        out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size + ei_idx*4+2] = s[1];
                        out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size + ei_idx*4+3] = s[0];
                    }
                    break;
                }
                case ECC3: {
                    memset(&out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size], 0xFF, oob_area_size);
                    for (ei_idx = 0 ; ei_idx < page_size/128 ; ei_idx++) {
                        unsigned char s[4] = {0};
                        reed_solomon_128bytes_ecc(&block_buffer[ p_idx*page_size + ei_idx*128], s, ECC3);
                        out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size + ei_idx*3  ] = s[2];
                        out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size + ei_idx*3+1] = s[1];
                        out_buffer[(p_idx+1)*page_size+p_idx*oob_area_size + ei_idx*3+2] = s[0];
                    }
                    break;
                }
            }
        }
        fwrite(out_buffer, block_size+(oob_area_size*pages_in_block), 1, out_file);
    }
    
exit:
    fclose(in_file);
    fclose(out_file);
    exit(0);
}