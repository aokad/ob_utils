import os, sys
import subprocess, shutil
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_scaffold
from .filter_bedpe import filter_scaffold_singleend_SV
from .repair_bedpe import repair_sv_type
from .svtools.vcftobedpe import run_vcf2bedpe 
from .svtools.vcftobed import run_vcf2bed 

def gridssSVtoBedpe(input_vcf, output, margin, f_grc, bcf_filter_option, debug, filter_scaffold_option):

    out_pref, ext = os.path.splitext(output)

    if bcf_filter_option != "":
        subprocess.check_call(['bcftools', 'view', '-i', 'INFO/MATEID[0]!=""', '-o', out_pref + '.tmp1.vcf', '-f', bcf_filter_option, input_vcf])
    else:
        subprocess.check_call(['bcftools', 'view', '-i', 'INFO/MATEID[0]!=""', '-o', out_pref + '.tmp1.vcf', input_vcf])
        
    
    run_vcf2bedpe(out_pref + ".tmp1.vcf", out_pref + ".tmp1.bedpe")

    give_margin_bedpe(out_pref + '.tmp1.bedpe', out_pref + '.tmp2.bedpe', margin)
    
    if filter_scaffold_option:
        filter_scaffold(out_pref + '.tmp2.bedpe', out_pref + '.tmp3.bedpe', f_grc)
    else:
        shutil.copyfile(out_pref + '.tmp2.bedpe', out_pref + '.tmp3.bedpe')
    
    hOUT = open(out_pref + '.tmp4.bedpe', 'w')
    subprocess.check_call(['bedtools', 'sort', '-i', out_pref + '.tmp3.bedpe'], stdout = hOUT)
    hOUT.close()

    repair_sv_type(out_pref + '.tmp4.bedpe', output)

    if not debug:
        os.remove(out_pref + '.tmp1.vcf')
        os.remove(out_pref + '.tmp1.bedpe')
        os.remove(out_pref + '.tmp2.bedpe')
        os.remove(out_pref + '.tmp3.bedpe')
        os.remove(out_pref + '.tmp4.bedpe')


def gridssSingleendSVtoBed(input_vcf, output, margin, f_grc, bcf_filter_option, debug, filter_scaffold_option):

    out_pref, ext = os.path.splitext(output)

    if bcf_filter_option != "":
        subprocess.check_call(['bcftools', 'view', '-e', 'INFO/MATEID[0]!=""', '-o', out_pref + '.tmp1.vcf', '-f', bcf_filter_option, input_vcf])
    else:
        subprocess.check_call(['bcftools', 'view', '-e', 'INFO/MATEID[0]!=""', '-o', out_pref + '.tmp1.vcf', input_vcf])
    
    run_vcf2bed(out_pref + ".tmp1.vcf", out_pref + ".tmp1.bed")

    if filter_scaffold_option:
        filter_scaffold_singleend_SV(out_pref + '.tmp1.bed', out_pref + '.tmp2.bed', f_grc)
    else:
        shutil.copyfile(out_pref + '.tmp1.bed', out_pref + '.tmp2.bed')

    hOUT = open(output, 'w')
    subprocess.check_call(['bedtools', 'sort', '-i', out_pref + '.tmp2.bed'], stdout = hOUT)
    hOUT.close()

    if not debug:
        os.remove(out_pref + '.tmp1.vcf')
        os.remove(out_pref + '.tmp1.bed')
        os.remove(out_pref + '.tmp2.bed')

def gridssSVtoBedpe_main(args):

    if args.singleend_sv:
        gridssSingleendSVtoBed(args.in_gridss_sv, args.output, args.margin, args.f_grc, args.bcf_filter_option, args.debug, args.filter_scaffold_option)
    else:
        gridssSVtoBedpe(args.in_gridss_sv, args.output, args.margin, args.f_grc, args.bcf_filter_option, args.debug, args.filter_scaffold_option)

