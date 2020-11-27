import os, sys
import subprocess, shutil
from . import utils
from .svtools.vcftobedpe import run_vcf2bedpe 
from .filter_bedpe import filter_scaffold
import pysam

def svimSVtoBedpe(input_vcf, output, f_grc, filter_scaffold_option, bcf_filter_option):

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

def filt_clustered_rearrangement2(input_file, output_file, control_junction_bedpe, control_check_margin, min_support_read,h_chrom_number):

    hout = open(output_file, 'w')
    control_junction_db = pysam.TabixFile(control_junction_bedpe)

    l_key = []    
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tstart1, tend1, tchr2, tstart2, tend2 = F[0], str(int(F[1])-1), F[2], F[3], str(int(F[4])-1), F[5]
            tdir1, tdir2 = F[8], F[9]
            infoA = F[18]
            support_read = utils.get_info_val(infoA, "SUPPORT")
            if int(support_read) < min_support_read: continue

            sort_flag = utils.sort_breakpoint_main(tchr1,tstart1,tchr2,tstart2,h_chrom_number)
            if not sort_flag:
                tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2 = tchr2, tstart2, tend2, tdir2, tchr1, tstart1, tend1, tdir1
            current_key = tchr1+"\t"+tend1+"\t"+tdir1+"\t"+tchr2+"\t"+tend2+"\t"+tdir2
        
            if current_key in l_key: continue

            control_flag = False
            tabix_error_flag = False
            try:
                records = control_junction_db.fetch(F[0], max(0, int(tstart1) - 200), int(tend1) + 200)
            except:
                tabix_error_flag = True

            if not tabix_error_flag:
                for record_line in records:
                    record = record_line.split('\t')

                    if tchr1 == record[0] and tdir1 == record[8] and \
                    int(tend1) >= int(record[1]) - control_check_margin and \
                    int(tstart1) <= int(record[2]) + control_check_margin and \
                    tchr2 == record[3] and tdir2 == record[9] and \
                    int(tend2) >= int(record[4]) - control_check_margin and \
                    int(tstart2) <= int(record[5]) + control_check_margin:
                        control_flag = True

            if not control_flag:
                print('\t'.join(F), file = hout)

            l_key.append(current_key)
    hout.close()
    
    

def simplify_svim(in_control_bedpe, min_support_read, hout, h_chrom_number):

    l_key = []
    with open(in_control_bedpe, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            tchr1, tstart1, tend1, tchr2, tstart2, tend2 = F[0], str(int(F[1])-1), F[2], F[3], str(int(F[4])-1), F[5]
            tdir1, tdir2 = F[8], F[9]
            info1 = F[18]
            support_read = utils.get_info_val(info1, "SUPPORT")
            if int(support_read) < min_support_read: continue

            sort_flag = utils.sort_breakpoint_main(tchr1,tstart1,tchr2,tstart2,h_chrom_number)
            if not sort_flag:
                tchr1, tstart1, tend1, tdir1, tchr2, tstart2, tend2, tdir2 = tchr2, tstart2, tend2, tdir2, tchr1, tstart1, tend1, tdir1

            current_key =tchr1+'\t'+tstart1+'\t'+tend1+'\t'+tdir1+'\t'+tchr2+'\t'+tstart2+'\t'+tend2+'\t'+tdir2

            if current_key in l_key:
                continue
            
            l_bed_record = [tchr1, str(int(tstart1)-1), tend1, tchr2, str(int(tstart2)-1), tend2, ".", ".", tdir1, tdir2]
            print('\t'.join(l_bed_record), file = hout)
            l_key.append(current_key)
            
            
def svimSVtoBedpe_main(args):
    
    in_tumor_sv = args.in_svim_tumor_sv
    in_control_sv = args.in_svim_control_sv
    margin = args.margin
    f_grc = args.f_grc
    filter_scaffold_option = args.filter_scaffold_option
    bcf_filter_option = args.bcf_filter_option
    min_control_support_read = args.min_control_support_read
    min_tumor_support_read = args.min_tumor_support_read

    output_prefix, ext = os.path.splitext(args.output)
    
    svimSVtoBedpe(in_tumor_sv, output_prefix+'.svim_tumor_PASS.bedpe', f_grc, filter_scaffold_option, bcf_filter_option)

    svimSVtoBedpe(in_control_sv, output_prefix+'.svim_control_PASS.bedpe', f_grc, filter_scaffold_option, bcf_filter_option)

    bcftools_command = ["bcftools", "view", "-h", in_tumor_sv, "-o", output_prefix +'.svim.vcf.header']
    subprocess.check_call(bcftools_command)
    h_chrom_number = utils.make_chrom_number_dict(output_prefix +'.svim.vcf.header')

    with open(output_prefix+'.svim_control_simplify.bedpe', 'w') as hout:
        simplify_svim(output_prefix+'.svim_control_PASS.bedpe', min_control_support_read, hout, h_chrom_number)

    with open( output_prefix +'.svim_control_sorted.bedpe', 'w') as hout:
        subprocess.check_call(['sort', '-k1,1', '-k2,2n', '-k4,4', '-k5,5n',  output_prefix +'.svim_control_simplify.bedpe'],  stdout = hout)

    with open(output_prefix +'.svim_control_sorted.bedpe.gz', "w") as hout:
        subprocess.check_call(["bgzip", "-f", "-c", output_prefix +'.svim_control_sorted.bedpe'], stdout = hout)
    subprocess.check_call(["tabix", "-p", "bed", output_prefix +'.svim_control_sorted.bedpe.gz'])
       
    filt_clustered_rearrangement2(output_prefix+'.svim_tumor_PASS.bedpe', output_prefix+'.svim_filtered.bedpe', 
    output_prefix+'.svim_control_sorted.bedpe.gz', margin, min_tumor_support_read, h_chrom_number)
    
    