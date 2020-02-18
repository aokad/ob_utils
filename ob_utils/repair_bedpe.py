from __future__ import print_function
import os, sys


def repair_dup_strand(bedpe_file, output):
    
    hOUT = open(output, 'w')
    with open(bedpe_file, 'r') as hin:
        for line in hin:
            
            if line.startswith("#"):
                header = line.rstrip('\n')
                print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')

            strandA = F[8]
            strandB = F[9]
            sv_type = F[10]
           
            if sv_type.startswith("DUP"):
                strandA = '-'
                strandB = '+'
            print('\t'.join(F[0:8])+'\t'+strandA+'\t'+strandB+'\t'+ '\t'.join(F[10:]), file=hOUT)
    hOUT.close()     
    

def repair_bedpe_main(args):

    repair_dup_strand(args.in_bedpe, args.output)
