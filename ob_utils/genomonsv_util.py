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


def genomonSVFormattoBedpe2(genomonsv_file, output):
   
    idx = 1
    d_svtype = {'deletion':'DEL', 'inversion':'INV', 'tandem_duplication':'DUP', 'translocation':'TRA'} 
    with open(output, 'w') as hout:
         with open(genomonsv_file, 'r') as hin:
             for line in hin:
                 if line.startswith("#"): continue
                 if line.startswith("Chr_1"):
                     header = line.rstrip('\n')
                     hF = header.split("\t")
                     print("\t".join(['#Chr1','Start1','End1','Chr2','Start2','End2', \
                           'Name','Score','Strand1','Strand2','Variant_type','Insert_seq','Genomonsv_info']), file=hout)
                     continue
                 line = line.rstrip('\n')
                 F = line.split('\t')
                 strand1, strand2 = F[2], F[5]
                 insseq = F[6]
                 svtype = F[7]
                 info = ";".join(['Num_Tumor_Ref_Read_Pair='+ F[12], \
                        'Num_Tumor_Var_Read_Pair='+ F[13], \
                        'Tumor_VAF='+ F[14], \
                        'Num_Control_Ref_Read_Pair='+ F[15], \
                        'Num_Control_Var_Read_Pair='+ F[16], \
                        'Control_VAF='+ F[17], \
                        'Minus_Log_Fisher_P_value='+ F[18], \
                        'Non-Matched_Control_Sample_With_Max_Junction='+ F[19], \
                        'Num_Max_Non-Matched_Control_Junction='+ F[20], \
                        'Max_Over_Hang_1='+ F[21], \
                        'Max_Over_Hang_2='+ F[22]])
                 print("\t".join([F[0],str(int(F[1])-1),F[1],F[3],str(int(F[4])-1),F[4], \
                       'GenomonSV_'+str(idx), 'NA', strand1, strand2, d_svtype[svtype], insseq, info]), file=hout)
                 idx +=1 

def genomonSVtoBedpe2(in_genomonsv, output, margin, f_grc):

    out_pref, ext = os.path.splitext(output)
    genomonSVFormattoBedpe2(in_genomonsv, out_pref + ".tmp1.bedpe")

    give_margin_bedpe(out_pref + ".tmp1.bedpe", out_pref + ".tmp2.bedpe", margin)
    filter_scaffold(out_pref + ".tmp2.bedpe", out_pref + ".tmp3.bedpe", f_grc)

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp3.bedpe"], stdout = hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.bedpe")
    os.remove(out_pref + ".tmp2.bedpe")
    os.remove(out_pref + ".tmp3.bedpe")

def genomonSVtoBedpe_main(args):
  
    if args.v2:  
        genomonSVtoBedpe2(args.in_genomon_sv, args.output, args.margin, args.f_grc)
    else:
        genomonSVtoBedpe(args.in_genomon_sv, args.output, args.margin, args.f_grc)

    
