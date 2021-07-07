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
    

def repair_sv_type(bedpe_file, output):
    
    hOUT = open(output, 'w')
    with open(bedpe_file, 'r') as hin:
        for line in hin:
            
            if line.startswith("#"):
                header = line.rstrip('\n')
                print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')

            chromA = F[0]
            chromB = F[3]
            strandA = F[8]
            strandB = F[9]
            sv_type = F[10]
           
            if chromA == chromB and strandA == '-' and strandB == '+':
                sv_type = "DUP"
            elif chromA == chromB and strandA == '+' and strandB == '-':
                sv_type = "DEL"
            elif chromA == chromB and strandA == '+' and strandB == '+':
                sv_type = "INV"
            elif chromA == chromB and strandA == '-' and strandB == '-':
                sv_type = "INV"

            print('\t'.join(F[0:10])+'\t'+sv_type+'\t'+ '\t'.join(F[11:]), file=hOUT)
    hOUT.close()    


def repair_dup_strand_cnv(bedpe_file, output):
    
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
            alt = F[14]
           
            if alt.startswith("<DUP>"):
                strandA = '-'
                strandB = '+'
            print('\t'.join(F[0:8])+'\t'+strandA+'\t'+strandB+'\t'+ '\t'.join(F[10:]), file=hOUT)
    hOUT.close()     

def repair_bedpe_main(args):

    repair_dup_strand(args.in_bedpe, args.output)
