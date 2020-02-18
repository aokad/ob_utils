from __future__ import print_function
import os, sys
import subprocess
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_scaffold
from .repair_bedpe import repair_dup_strand

def SvABASVtoBedpe(input_vcf, output, margin, f_grc, normal_max_variant_read, 
    tumor_min_variant_read, tumor_min_depth, normal_min_depth, f_germ):

    out_pref, ext = os.path.splitext(output)

    subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", "-f", "PASS", input_vcf])

    command = ['bcftools', 'filter', '--output', out_pref + '.tmp2.vcf', '--include']
    if f_germ:
        command.append('FORMAT/AD[0:0]>='+str(tumor_min_variant_read)+'&&FORMAT/DP[0:0]>='+str(tumor_min_depth))
    else:
        command.append('FORMAT/AD[0:0]<='+str(normal_max_variant_read)+'&&FORMAT/AD[1:0]>='+str(tumor_min_variant_read)+'&&FORMAT/DP[0:0]>='+str(normal_min_depth)+'&&FORMAT/DP[1:0]>='+str(tumor_min_depth))
    command.append(out_pref + ".tmp1.vcf")
    subprocess.check_call(command)

    subprocess.check_call(["svtools", "vcftobedpe", "-i", out_pref + ".tmp2.vcf", "-o", out_pref + ".tmp2.bedpe"])

    give_margin_bedpe(out_pref + ".tmp2.bedpe", out_pref + ".tmp3.bedpe", margin)
    filter_scaffold(out_pref + ".tmp3.bedpe", out_pref + ".sv.bedpe", f_grc)

    os.remove(out_pref + ".tmp1.vcf")
    os.remove(out_pref + ".tmp2.vcf")
    os.remove(out_pref + ".tmp2.bedpe")
    os.remove(out_pref + ".tmp3.bedpe")


def filter_indel_vcf(in_vcf, outfile, min_del_size, min_ins_size):

    f_header = True
    hOUT = open(outfile, 'w')
    with open(in_vcf, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            
            if line.startswith("##"): # if line is the meta-data
                print(line, file=hOUT)
                continue
            elif line.startswith("#"): # if line is the header-line
                if f_header: # add the original metadata
                    print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">', file=hOUT)
                    print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">', file=hOUT)
                    f_header = False
                print(line, file=hOUT)
                continue
            
            F = line.split('\t')
            POS = int(F[1])
            REF = F[3]
            ALT = F[4]
            INFO = F[7] 

            SVTYPE = 'DUP|INS' if len(REF) < len(ALT) else 'DEL'

            END = POS
            if SVTYPE == 'DUP|INS':
                if (len(ALT) - len(REF)) < min_ins_size: continue
                END += (len(ALT) - len(REF))
            elif SVTYPE == 'DEL':
                if (len(REF) - len(ALT)) < min_del_size: continue
                END += (len(REF) - len(ALT))

            new_info = INFO+';SVTYPE='+SVTYPE+';END='+str(END)
            print('\t'.join(F[:7]) + '\t'+new_info+'\t'+ '\t'.join(F[8:]), file=hOUT)
        

def SvABAIndeltoBedpe(input_vcf, output, margin, f_grc, normal_max_variant_read,
    tumor_min_variant_read, tumor_min_depth, normal_min_depth, min_del_size, min_ins_size, f_germ):

    out_pref, ext = os.path.splitext(output)

    subprocess.check_call(["bcftools", "view", "-o", out_pref + ".indel.tmp1.vcf", "-f", "PASS", input_vcf])

    command = ['bcftools', 'filter', '--output', out_pref + '.indel.tmp2.vcf', '--include']
    if f_germ:
        command.append('FORMAT/AD[0:0]>='+str(tumor_min_variant_read)+'&&FORMAT/DP[0:0]>='+str(tumor_min_depth))
    else:
        command.append('FORMAT/AD[0:0]<='+str(normal_max_variant_read)+'&&FORMAT/AD[1:0]>='+str(tumor_min_variant_read)+'&&FORMAT/DP[0:0]>='+str(normal_min_depth)+'&&FORMAT/DP[1:0]>='+str(tumor_min_depth))
    command.append(out_pref + ".indel.tmp1.vcf")
    subprocess.check_call(command)
    
    filter_indel_vcf(out_pref + ".indel.tmp2.vcf", out_pref + ".indel.tmp3.vcf", min_del_size, min_ins_size)

    subprocess.check_call(["svtools", "vcftobedpe", "-i", out_pref + ".indel.tmp3.vcf", "-o", out_pref + ".indel.tmp3.bedpe"])

    give_margin_bedpe(out_pref + ".indel.tmp3.bedpe", out_pref + ".indel.tmp4.bedpe", margin)
    filter_scaffold(out_pref + ".indel.tmp4.bedpe", out_pref + ".indel.tmp5.bedpe", f_grc)
    repair_dup_strand(out_pref + ".indel.tmp5.bedpe", out_pref + ".indel.bedpe")


    os.remove(out_pref + ".indel.tmp1.vcf")
    os.remove(out_pref + ".indel.tmp2.vcf")
    os.remove(out_pref + ".indel.tmp3.vcf")
    os.remove(out_pref + ".indel.tmp3.bedpe")
    os.remove(out_pref + ".indel.tmp4.bedpe")
    os.remove(out_pref + ".indel.tmp5.bedpe")


def merge_sv_and_indel_bedpe(in_sv_bedpe,in_indel_bedpe, output):

    out_pref, ext = os.path.splitext(output)
    
    hOUT = open(out_pref+".tmp1.bedpe", 'w')
    with open(in_sv_bedpe, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n') 
            print(line, file=hOUT)

    with open(in_indel_bedpe, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n') 
            if line.startswith('#'): continue
            print(line, file=hOUT)
    hOUT.close()

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref+".tmp1.bedpe"], stdout=hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.bedpe")
    os.remove(in_sv_bedpe)
    os.remove(in_indel_bedpe)


def svabaSVtoBedpe_main(args):

    SvABASVtoBedpe(args.in_svaba_sv, args.output, args.margin, args.f_grc, 
    args.normal_max_variant, args.tumor_min_variant, 
    args.tumor_min_depth, args.normal_min_depth, args.f_germ)

    SvABAIndeltoBedpe(args.in_svaba_indel, args.output, args.margin, args.f_grc, 
    args.normal_max_variant, args.tumor_min_variant, 
    args.tumor_min_depth, args.normal_min_depth,
    args.min_del_size, args.min_ins_size, args.f_germ)
    
    out_pref, ext = os.path.splitext(args.output)
    merge_sv_and_indel_bedpe(out_pref + ".sv.bedpe", out_pref + ".indel.bedpe", args.output)
    
    