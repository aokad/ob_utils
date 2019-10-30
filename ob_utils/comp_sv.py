from __future__ import print_function
import os, sys
import subprocess


def genomonSVtoBedpe(genomonsv_file, margin, output):
    
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
            startA = str(int(F[1]) - margin)
            endA = str(int(F[1]) + margin)
            strandA = F[2]
            posA = F[1]
            chrB = F[3]
            startB = str(int(F[4]) - margin)
            endB = str(int(F[4]) + margin)
            strandB = F[5]
            posB = F[4]
            print(chrA+'\t'+startA+'\t'+endA+'\t'+chrB+'\t'+startB+'\t'+endB+'\t.\t.\t'+strandA+'\t'+strandB+'\t'+posA+'\t'+posB+'\t'+ '\t'.join(F[6:]), file=hOUT)
    hOUT.close()


def onebreaktoBed(onebreak_file, margin, output):
    
    hOUT = open(output, 'w')
    with open(onebreak_file, 'r') as hin:
        for line in hin:
            if line.startswith("Chr"):
                header = line.rstrip('\n')
                hF = header.split("\t")
                print('#CHROM\tSTART\tEND\tID\tQUAL\tSTRAND\tPOS\t'+ '\t'.join(hF[3:]), file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrA = F[0]
            startA = str(int(F[1]) - margin)
            endA = str(int(F[1]) + margin)
            strandA = F[2]
            posA = F[1]
            print(chrA+'\t'+startA+'\t'+endA+'\t.\t.\t'+strandA+'\t'+posA+'\t'+ '\t'.join(F[3:]), file=hOUT)
        hOUT.close()


def makeHash(anno_file, idx1, idx2, idx3, margin):
    sv_comp = {}
    with open(anno_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            chrA = F[0]
            chrB = F[3]
            strandA = F[8]
            strandB = F[9]
            posA = str(int(F[1]) + margin)
            posB = str(int(F[4]) + margin)
            chrC = F[idx1]
            posC = F[idx2]
            strandC = F[idx3]

            key = chrC+'\t'+posC+'\t'+strandC
            val = chrA+':'+strandA+posA+'-'+chrB+':'+strandB+posB
            sv_comp[key] = val
    return sv_comp


def annotSvComp(onebreak_file, sv_comp1, output):
    
    hOUT = open(output, 'w')
    with open(onebreak_file, 'r') as hin:
        for line in hin:
            if line.startswith("Chr"):
                header = line.rstrip('\n')
                print(header+'\tGenomonSV', file=hOUT)
                continue
            F = line.rstrip('\n').split('\t')
            chrA = F[0]
            posA = F[1]
            strandA = F[2]
            key = chrA+'\t'+posA+'\t'+strandA
            
            sv_info1 = sv_comp1[key] if key in sv_comp1 else "---"
            print('\t'.join(F) +'\t'+ sv_info1, file=hOUT)
    hOUT.close()


def comp_main(args):
    
    out_pref, ext = os.path.splitext(args.output)
    genomonSVtoBedpe(args.in_genomonsv, args.margin, out_pref + ".sv1.bedpe")
    
    # hOUT = open(out_pref + ".sv2.bedpe", 'w')
    # subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".sv1.bedpe"], stdout = hOUT)
    # hOUT.close()
    
    onebreaktoBed(args.in_onebreak, args.margin, out_pref + ".ob1.bed")

    hOUT = open(out_pref + ".pair.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtobed", "-s", "-a", out_pref + ".sv1.bedpe", "-b", out_pref + ".ob1.bed"], stdout = hOUT)
    hOUT.close()
    
    sv_comp1 = makeHash(out_pref + ".pair.bedpe", 29, 35, 34, args.margin)
    
    annotSvComp(args.in_onebreak, sv_comp1, args.output)
    
    