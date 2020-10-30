from __future__ import print_function
import os, sys

def give_margin_bedpe(bedpe_file, output, margin):

    hOUT = open(output, 'w')
    with open(bedpe_file, 'r') as hin:
        for line in hin:
            if line.startswith("#"):
                header = line.rstrip('\n')
                print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrA = F[0]
            startA = str(max(0,int(F[1]) - margin))
            endA = str(int(F[2]) + margin)
            chrB = F[3]
            startB = str(max(0,int(F[4]) - margin))
            endB = str(int(F[5]) + margin)
            print(chrA+'\t'+startA+'\t'+endA+'\t'+chrB+'\t'+startB+'\t'+endB+'\t'+ '\t'.join(F[6:]), file=hOUT)
    hOUT.close()       

def margin_bedpe_main(args):

    give_margin_bedpe(args.in_bedpe, args.output, args.margin)
    