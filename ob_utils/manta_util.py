import os, sys
import subprocess
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_scaffold
from .repair_bedpe import repair_dup_strand

def mantaSVtoBedpe(input_vcf, output, margin, f_grc):

    out_pref, ext = os.path.splitext(output)

    subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", "-f", "PASS", input_vcf])

    subprocess.check_call(["svtools", "vcftobedpe", "-i", out_pref + ".tmp1.vcf", "-o", out_pref + ".tmp1.bedpe"])

    give_margin_bedpe(out_pref + ".tmp1.bedpe", out_pref + ".tmp2.bedpe", margin)
    filter_scaffold(out_pref + ".tmp2.bedpe", out_pref + ".tmp3.bedpe", f_grc)
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

    mantaSVtoBedpe(args.in_manta_sv, args.output, args.margin, args.f_grc)

    
