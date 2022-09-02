import os, sys
import subprocess, shutil
from . import utils
from .svtools.vcftobedpe import run_vcf2bedpe 
import pysam


def simplify_sniffles(in_control_bedpe, min_support_read, hout):

    with open(in_control_bedpe, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            chrA, chrB, strandA, strandB = F[0], F[3], F[8], F[9]
            infoA = F[18]
            posA, posB = utils.get_position(infoA, "")
            sv_type = F[10]
            alt = F[14]
            support_read = utils.get_info_val(infoA, "RE")
            if int(support_read) < min_support_read: continue

            l_bed_record = [chrA, str(int(posA)-1), posA, chrB, str(int(posB)-1), posB, ".", ".", strandA, strandB]
            print('\t'.join(l_bed_record), file = hout)


def snifflesSVtoBedpe_main(args):
    print("called sniffles_control:snifflesSVtoBedpe_main")
    withopen
    simplify_sniffles(args.in_bedpe, args.output, args.in_control_bedpe)

