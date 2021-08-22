from __future__ import print_function
import os, sys
import subprocess
from .margin_bedpe import give_margin_bedpe
from .filter_bedpe import filter_scaffold

def bedpetoBedpeForComparison(in_bedpe, output, margin):

    out_pref, ext = os.path.splitext(output)

    give_margin_bedpe(in_bedpe, out_pref + ".tmp1.bedpe", margin)

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp1.bedpe"], stdout = hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.bedpe")

def bedpetoBedpe_main(args):
  
    bedpetoBedpeForComparison(args.in_bedpe, args.output, args.margin)

    
