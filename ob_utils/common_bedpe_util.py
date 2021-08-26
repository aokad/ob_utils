from __future__ import print_function
import os, sys
import subprocess
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_min_sv_size

def bedpetoBedpeForComparison(in_bedpe, output, margin, min_sv_size):

    out_pref, ext = os.path.splitext(output)

    filter_min_sv_size(in_bedpe, out_pref+".tmp1.bedpe", min_sv_size)

    give_margin_bedpe(out_pref+".tmp1.bedpe", out_pref+".tmp2.bedpe", margin)

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp2.bedpe"], stdout = hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.bedpe")
    os.remove(out_pref + ".tmp2.bedpe")

def bedpetoBedpe_main(args):
  
    bedpetoBedpeForComparison(args.in_bedpe, args.output, args.margin, args.min_sv_size)

    
