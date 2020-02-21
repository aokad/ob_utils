import os, sys
import subprocess
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_scaffold
from .repair_bedpe import repair_sv_type

def gridssSVtoBedpe(input_vcf, output, margin, f_grc):

    out_pref, ext = os.path.splitext(output)

    subprocess.check_call(['bcftools', 'view', '-i', 'INFO/MATEID[0]!=""', '-o', out_pref + '.tmp1.vcf', '-f', 'PASS', input_vcf])
    
    # subprocess.check_call(['sed','-i',"s/PARID=/MATEID=/g", out_pref + '.tmp1.vcf'])

    subprocess.check_call(['svtools', 'vcftobedpe', '-i', out_pref + '.tmp1.vcf', '-o', out_pref + '.tmp1.bedpe'])

    give_margin_bedpe(out_pref + '.tmp1.bedpe', out_pref + '.tmp2.bedpe', margin)
    filter_scaffold(out_pref + '.tmp2.bedpe', out_pref + '.tmp3.bedpe', f_grc)
    
    hOUT = open(out_pref + '.tmp4.bedpe', 'w')
    subprocess.check_call(['bedtools', 'sort', '-i', out_pref + '.tmp3.bedpe'], stdout = hOUT)
    hOUT.close()

    repair_sv_type(out_pref + '.tmp4.bedpe', output)

    os.remove(out_pref + '.tmp1.vcf')
    os.remove(out_pref + '.tmp1.bedpe')
    os.remove(out_pref + '.tmp2.bedpe')
    os.remove(out_pref + '.tmp3.bedpe')
    os.remove(out_pref + '.tmp4.bedpe')

def gridssSVtoBedpe_main(args):

    gridssSVtoBedpe(args.in_gridss_sv, args.output, args.margin, args.f_grc)

    
