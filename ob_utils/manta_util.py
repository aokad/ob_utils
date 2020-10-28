import os, sys
import subprocess, shutil
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_scaffold
from .repair_bedpe import repair_dup_strand
from .svtools.vcftobedpe import run_vcf2bedpe 

def mantaSVtoBedpe(input_vcf, output, margin, f_grc, bcf_filter_option, filter_scaffold_option):

    out_pref, ext = os.path.splitext(output)

    if bcf_filter_option != "":
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", "-f", bcf_filter_option, input_vcf])
    else:
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", input_vcf])
    run_vcf2bedpe(out_pref + ".tmp1.vcf", out_pref + ".tmp1.bedpe")

    give_margin_bedpe(out_pref + ".tmp1.bedpe", out_pref + ".tmp2.bedpe", margin)
    
    if filter_scaffold_option:
        filter_scaffold(out_pref + ".tmp2.bedpe", out_pref + ".tmp3.bedpe", f_grc)
    else:
        shutil.copyfile(out_pref + '.tmp2.bedpe', out_pref + '.tmp3.bedpe')
        
    repair_dup_strand(out_pref + ".tmp3.bedpe", out_pref + ".tmp4.bedpe")

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp4.bedpe"], stdout = hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.vcf")
    os.remove(out_pref + ".tmp1.bedpe")
    os.remove(out_pref + ".tmp2.bedpe")
    os.remove(out_pref + ".tmp3.bedpe")
    os.remove(out_pref + ".tmp4.bedpe")

def mantaSVtoBedpe_main(args):
    mantaSVtoBedpe(args.in_manta_sv, args.output, args.margin, args.f_grc, args.bcf_filter_option, args.filter_scaffold_option)

    
