import os, sys
import subprocess, shutil
from . import utils
from .svtools.vcftobedpe import run_vcf2bedpe 
from .filter_bedpe import filter_scaffold
import pysam
import re

def snifflesSVtoBedpe(input_vcf, output, f_grc, filter_scaffold_option, bcf_filter_option):

    out_pref, ext = os.path.splitext(output)

    if bcf_filter_option != "":
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", "-f", bcf_filter_option, input_vcf])
    else:
        subprocess.check_call(["bcftools", "view", "-o", out_pref + ".tmp1.vcf", input_vcf])
        
    run_vcf2bedpe(out_pref + ".tmp1.vcf", out_pref + ".tmp1.bedpe")
    
    if filter_scaffold_option:
        filter_scaffold(out_pref + ".tmp1.bedpe", out_pref + ".tmp2.bedpe", f_grc)
    else:
        shutil.copyfile(out_pref + '.tmp1.bedpe', out_pref + '.tmp2.bedpe')

    hOUT = open(output, 'w')
    subprocess.check_call(["bedtools", "sort", "-i", out_pref + ".tmp2.bedpe"], stdout = hOUT)
    hOUT.close()

    os.remove(out_pref + ".tmp1.vcf")
    os.remove(out_pref + ".tmp1.bedpe")
    os.remove(out_pref + ".tmp2.bedpe")


def repair_dup_strand(bedpe_file, output):
    
    hOUT = open(output, 'w')
    with open(bedpe_file, 'r') as hin:
        for line in hin:
            
            if line.startswith("#"):
                header = line.rstrip('\n')
                print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')

            sv_type = F[10]
            alt = F[14]

            if sv_type == "DEL":
                F[8] = '+'
                F[9] = '-'
            elif sv_type == "INV":
                F[8] = '*'
                F[9] = '*'
            elif sv_type == "DUP":
                F[8] = '-'
                F[9] = '+'
            elif sv_type == "INS":
                F[8] = '+'
                F[9] = '-'
            #elif sv_type == "BND":
            #    if alt.startswith('N'):
            #        F[8] = '+'
            #    else:
            #        F[8] = '-'
            #    alt_seq = re.findall(r'([][])(.+?)([][])', alt)
            #    if alt_seq[0][0].startswith(']'):
            #        F[9] = '+'
            #    else:
            #        F[9] = '-'
            print('\t'.join(F), file=hOUT)
            
    hOUT.close()  

def filt_clustered_rearrangement2(input_file, output_file, control_junction_bedpe, control_check_margin, min_tumor_support_read,max_control_support_read,min_sv_length,h_chrom_number, sniffles2):

    hout = open(output_file, 'w')
    control_junction_db = pysam.TabixFile(control_junction_bedpe)
    
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tstart1, tend1, tchr2, tstart2, tend2 = F[0], str(int(F[1])-1), F[2], F[3], str(int(F[4])-1), F[5]
            tdir1, tdir2 = F[8], F[9]
            sv_type = F[10]
            alt = F[14]
            info1, info2 = F[18], F[19]
            format_keys, format_vals = F[20], F[21]
            if sniffles2:
                tumor_support_read = utils.get_info_val(info1, "SUPPORT")
            else:
                tumor_support_read = utils.get_info_val(info1, "RE")
            if int(tumor_support_read) < min_tumor_support_read: continue
            # if sv_type != "BND":
            #    sv_len = abs(int(utils.get_info_val(info1, "SVLEN")))
            #    if sv_len < min_sv_length: continue        
            tumor_ref_read = utils.get_format_val(format_keys, format_vals, "DR")
            insert_seq = alt if sv_type == "INS" else "---"
                            
            sort_flag = utils.sort_breakpoint_main(tchr1,tstart1,tchr2,tstart2,h_chrom_number)
            if not sort_flag:
                tchr1, tstart1, tend1, tdir1, info1, tchr2, tstart2, tend2, tdir2, info2 = tchr2, tstart2, tend2, tdir2, info2, tchr1, tstart1, tend1, tdir1, info1

            control_flag = False
            tabix_error_flag = False
            try:
                records = control_junction_db.fetch(F[0], max(0, int(tstart1) - 200), int(tend1) + 200)
            except:
                tabix_error_flag = True

            control_ref_read = 0
            control_support_read = 0
            if not tabix_error_flag:
                for record_line in records:
                    record = record_line.split('\t')

                    if tchr1 == record[0] and tdir1 == record[8] and \
                    int(tstart1) >= int(record[1]) - control_check_margin and \
                    int(tstart1) <= int(record[1]) + control_check_margin and \
                    int(tend1) >= int(record[2]) - control_check_margin and \
                    int(tend1) <= int(record[2]) + control_check_margin and \
                    tchr2 == record[3] and tdir2 == record[9] and \
                    int(tstart2) >= int(record[4]) - control_check_margin and \
                    int(tstart2) <= int(record[4]) + control_check_margin and \
                    int(tend2) >= int(record[5]) - control_check_margin and \
                    int(tend2) <= int(record[5]) + control_check_margin:

                        control_support_read = int(record[11]) if record[11] != "." else record[11]
                        control_ref_read = int(record[12]) if record[12] != "." else record[12]
                        
                        if int(record[11]) > max_control_support_read:
                            control_flag = True
                            break

            if not control_flag:
                print("\t".join([tchr1, tend1, tdir1, tchr2, tend2, tdir2, insert_seq, tumor_ref_read, tumor_support_read, str(control_ref_read), str(control_support_read),sv_type]), file = hout)

    hout.close()

    
def simplify_sniffles(in_control_bedpe, hout, h_chrom_number, sniffles2):

    with open(in_control_bedpe, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            tchr1, tstart1, tend1, tchr2, tstart2, tend2, tdir1, tdir2 = F[0], F[1], F[2], F[3], F[4], F[5], F[8], F[9]
            sv_type = F[10]
            info1 = F[18]
            format_keys, format_vals = F[20], F[21]
            if sniffles2:
                support_read = utils.get_info_val(info1, "SUPPORT")
            else:
                support_read = utils.get_info_val(info1, "RE")
            ref_read = utils.get_format_val(format_keys, format_vals, "DR")
        
            sort_flag = utils.sort_breakpoint_main(tchr1,tstart1,tchr2,tstart2,h_chrom_number)
            if not sort_flag:
                tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2 = tchr2, tstart2, tend2, tdir2, tchr1, tstart1, tend1, tdir1

            l_bed_record = [tchr1, str(int(tstart1)-1), tend1, tchr2, str(int(tstart2)-1), tend2, ".", ".", tdir1, tdir2, sv_type, support_read, ref_read]
            print('\t'.join(l_bed_record), file = hout)
            
            
def snifflesSVtoBedpe_main(args):
    
    in_tumor_sv = args.in_sniffles_tumor_sv
    in_control_sv = args.in_sniffles_control_sv
    margin = args.margin
    f_grc = args.f_grc
    filter_scaffold_option = args.filter_scaffold_option
    bcf_filter_option = args.bcf_filter_option
    max_control_support_read = args.max_control_support_read
    min_tumor_support_read = args.min_tumor_support_read
    min_sv_length = args.min_sv_length
    output = args.output
    debug = args.debug
    sniffles2 = args.sniffles2

    output_prefix, ext = os.path.splitext(output)
    
    snifflesSVtoBedpe(in_tumor_sv, output_prefix+'.sniffles_tumor_PASS.bedpe', f_grc, filter_scaffold_option, bcf_filter_option)

    repair_dup_strand(output_prefix+'.sniffles_tumor_PASS.bedpe', output_prefix+'.sniffles_tumor_repaired.bedpe')

    snifflesSVtoBedpe(in_control_sv, output_prefix+'.sniffles_control_PASS.bedpe', f_grc, filter_scaffold_option, bcf_filter_option)

    repair_dup_strand(output_prefix+'.sniffles_control_PASS.bedpe', output_prefix+'.sniffles_control_repaired.bedpe')

    bcftools_command = ["bcftools", "view", "-h", in_tumor_sv, "-o", output_prefix +'.sniffles.vcf.header']
    subprocess.check_call(bcftools_command)
    h_chrom_number = utils.make_chrom_number_dict(output_prefix +'.sniffles.vcf.header')

    with open(output_prefix+'.sniffles_control_simplify.bedpe', 'w') as hout:
        simplify_sniffles(output_prefix+'.sniffles_control_repaired.bedpe', hout, h_chrom_number, sniffles2)

    with open( output_prefix +'.sniffles_control_sorted.bedpe', 'w') as hout:
        subprocess.check_call(['sort', '-k1,1', '-k2,2n', '-k4,4', '-k5,5n', '-k9,9', '-k10,10',  output_prefix +'.sniffles_control_simplify.bedpe'],  stdout = hout)

    with open(output_prefix +'.sniffles_control_sorted.bedpe.gz', "w") as hout:
        subprocess.check_call(["bgzip", "-f", "-c", output_prefix +'.sniffles_control_sorted.bedpe'], stdout = hout)
    subprocess.check_call(["tabix", "-p", "bed", output_prefix +'.sniffles_control_sorted.bedpe.gz'])
          
    filt_clustered_rearrangement2(output_prefix+'.sniffles_tumor_repaired.bedpe', output_prefix+'.sniffles_filtered.txt', 
    output_prefix+'.sniffles_control_sorted.bedpe.gz', margin, min_tumor_support_read, max_control_support_read, min_sv_length, h_chrom_number, sniffles2)

    with open(output_prefix + ".sniffles_sorted.txt", 'w') as hout:
        subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", "-k3,3", "-k6,6", output_prefix + ".sniffles_filtered.txt"],  stdout = hout)

    with open(output, 'w') as hout:
        l_header = ["Chr_1","Pos_1","Dir_1","Chr_2","Pos_2","Dir_2","Inserted_Seq","Checked_Read_Num_Tumor","Supporting_Read_Num_Tumor","Checked_Read_Num_Control","Supporting_Read_Num_Control","Sv_Type"]
        print("\t".join(l_header), file=hout)
        with open(output_prefix + ".sniffles_sorted.txt", 'r') as hin:
            for line in hin:
                print(line.rstrip('\n'), file=hout)
            
    if not debug:
        os.remove(output_prefix +'.sniffles_tumor_PASS.bedpe')
        os.remove(output_prefix +'.sniffles_control_PASS.bedpe')
        os.remove(output_prefix +'.sniffles_tumor_repaired.bedpe')
        os.remove(output_prefix +'.sniffles_control_repaired.bedpe')
        os.remove(output_prefix +'.sniffles.vcf.header')
        os.remove(output_prefix +'.sniffles_control_simplify.bedpe')
        os.remove(output_prefix +'.sniffles_control_sorted.bedpe')
        os.remove(output_prefix +'.sniffles_filtered.txt')
        os.remove(output_prefix +'.sniffles_sorted.txt')

    
