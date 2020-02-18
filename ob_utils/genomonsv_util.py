from __future__ import print_function
import os, sys
import subprocess
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_scaffold

def genomonSVFormattoBedpe(genomonsv_file, output):
    
    hOUT = open(output, 'w')
    with open(genomonsv_file, 'r') as hin:
        for line in hin:
            if line.startswith("#"): continue
            if line.startswith("Chr_1"):
                header = line.rstrip('\n')
                hF = header.split("\t")
                print('#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tPOS_A\tPOS_B\t'+ '\t'.join(hF[6:]), file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrA = F[0]
            startA = str(int(F[1]))
            endA = str(int(F[1]))
            strandA = F[2]
            posA = F[1]
            chrB = F[3]
            startB = str(int(F[4]))
            endB = str(int(F[4]))
            strandB = F[5]
            posB = F[4]
            print(chrA+'\t'+startA+'\t'+endA+'\t'+chrB+'\t'+startB+'\t'+endB+'\t.\t.\t'+strandA+'\t'+strandB+'\t'+posA+'\t'+posB+'\t'+ '\t'.join(F[6:]), file=hOUT)
    hOUT.close()
    

def genomonSVtoBedpe(in_genomonsv, output, margin, f_grc):

    out_pref, ext = os.path.splitext(output)
    genomonSVFormattoBedpe(in_genomonsv, out_pref + ".tmp1.bedpe")

    give_margin_bedpe(out_pref + ".tmp1.bedpe", out_pref + ".tmp2.bedpe", margin)
    filter_scaffold(out_pref + ".tmp2.bedpe", out_pref + ".tmp3.bedpe", f_grc)

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp3.bedpe"], stdout = hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.bedpe")
    os.remove(out_pref + ".tmp2.bedpe")
    os.remove(out_pref + ".tmp3.bedpe")

def genomonSVtoBedpe_main(args):
    
    genomonSVtoBedpe(args.in_genomon_sv, args.output, args.margin, args.f_grc)

    
