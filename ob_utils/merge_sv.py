from __future__ import print_function 
import os, sys

def merge_main(args):

    d_file2_header = ""
    d_file2 = {}

    with open(args.in_onebreak_filt2, 'r') as hin:
        for line in hin:
            if line.startswith("#"):
                header = line.rstrip('\n')
                F = header.split('\t')
                d_file2_header = "\t".join(F[7:])
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrom = F[0]
            pos = F[1]
            dir = F[2]
            key = chrom +"\t"+ pos +"\t"+ dir
            d_file2[key] = "\t".join(F[7:])
    
    with open(args.output, 'w') as hout:
        with open(args.in_onebreak_filt1, 'r') as hin:
            for line in hin:
                if line.startswith("#"):
                    header = line.rstrip('\n')
                    print(header+"\t"+d_file2_header, file=hout)
                    continue
                line = line.rstrip('\n')
                F = line.split('\t')
                chrom = F[0]
                pos = F[1]
                dir = F[2]
                key = chrom +"\t"+ pos +"\t"+ dir
                if key in d_file2:
                    print(line +"\t"+ d_file2[key], file=hout)
        
