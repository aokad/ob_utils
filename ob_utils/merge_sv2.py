from __future__ import print_function 
import os, sys
import subprocess
import pysam


def get_inseq_from_bp(bp_pair):
    chr1, pos1, dir1, chr2, pos2, dir2, inseq = get_position_from_bp(bp_pair)
    return inseq 

def get_position_from_bp(bp_pair):
    # format 
    # chr1,pos1,dir1,chr2,pos2,dir2,inseq
    # 1,4669434,-,3,135016646,GG
    
    chr1, pos1, dir1, chr2, pos2, dir2, inseq = bp_pair.split(",")

    # chr1, pos1dir1 = bp1.split(":")
    # dir1 = pos1dir1[-1]
    # pos1 = pos1dir1[:-1]

    # chr2, pos2dir2inseq = bp2.split(":")
    # idx = pos2dir2inseq.find('+')
    # dir2 = '+' if idx != -1 else '-'
    # pos2, inseq = pos2dir2inseq.split(dir2)
    
    return chr1, pos1, dir1, chr2, pos2, dir2, inseq 


def check_bp(tool1, tool2):

    t1chr1, t1pos1, t1dir1, t1chr2, t1pos2, t1dir2, inseq1 = get_position_from_bp(tool1)
    t2chr1, t2pos1, t2dir1, t2chr2, t2pos2, t2dir2, inseq2 = get_position_from_bp(tool2)
            
    check_margin_size = 20
    flag = 0
    match = 0
    # detailed check on the junction position considering inserted sequences
    
    if t1chr1 != t2chr1 or int(t2pos1) > int(t1pos1) + check_margin_size:
        return match
    else:
        if t1chr1 == t2chr1 and t1chr2 == t2chr2 and t1dir1 == t2dir1 and t1dir2 == t2dir2:
                
            if t2dir1 == "+":
                expectedDiffSize = (int(t2pos1) - int(t1pos1)) + (int(len(inseq2)) - int(len(inseq1)))
                if (t2dir2 == "+" and int(t2pos2) == int(t1pos2) - int(expectedDiffSize)) or (t2dir2 == "-" and int(t2pos2) == int(t1pos2) + int(expectedDiffSize)):
                    flag = 1
                            
            else:
                expectedDiffSize = (int(t2pos1) - int(t1pos1)) + (int(len(inseq1)) - int(len(inseq2)))
                if (t2dir2 == "+" and int(t2pos2) == int(t1pos2) + int(expectedDiffSize)) or (t2dir2 == "-" and int(t2pos2) == int(t1pos2) - int(expectedDiffSize)):
                    flag = 1
                
            # if the junction position and direciton match
            if flag == 1:
                match = 1
            
    return match


def get_info_val(infos, key):
    ret = ""
    l_info = infos.split(';')
    for info in l_info:
        if info.startswith(key):
            ret = info.split("=")[1]
            break
    return ret

def get_alt_seq(alt):
    ret = ""
    if alt.find("[") > -1:
        l_alt = alt.split("[")
        if len(l_alt[0]) > 1:   ret = (l_alt[0])[1:]
        elif len(l_alt[2]) > 1: ret = (l_alt[2])[:-1]
            
    elif alt.find("]") > -1:
        l_alt = alt.split("]")
        if len(l_alt[0]) > 1:   ret = (l_alt[0])[1:]
        elif len(l_alt[2]) > 1: ret = (l_alt[2])[:-1]
        
    return ret
    
def get_insert_seq(file_type, F, column_idx, alt_idx):
    insert_seq = ""
    if file_type == "GenomonSV":
        insert_seq = F[column_idx] 
    elif file_type == "Manta":
        insert_seq = get_info_val(F[column_idx], "SVINSSEQ")
    elif file_type == "SvABA":
        insert_seq = get_info_val(F[column_idx], "INSERTION")
    elif file_type == "GRIDSS":
        insert_seq = get_alt_seq(F[alt_idx])
    if insert_seq == "---": insert_seq = ""
    return insert_seq
    
def get_sv_type(sv_type):
    
    if sv_type == 'deletion': sv_type = 'DEL'
    if sv_type == 'inversion': sv_type = 'INV'
    if sv_type == 'tandem_duplication': sv_type = 'DUP'
    if sv_type == 'translocation': sv_type = 'TRANS'
    return sv_type

def get_sv_type_from_BND(chr1, dir1, chr2, dir2):

    sv_type = "TRANS"
    if chr1 == chr2 and dir1 == dir2:
        sv_type = "INV"
    return sv_type

def getpos(info):
    position = ''
    l_info = info.split(';')
    for v in l_info:
        if v.startswith('POS='):
            position = v.replace('POS=','')
            break
    return position


def getposend(info):
    position = ''
    endpos = ''
    l_info = info.split(';')
    for v in l_info:
        if v.startswith('POS='):
            position = v.replace('POS=','')
        if v.startswith('END='):
            endpos = v.replace('END=','')
    return position, endpos
    

def get_position(infoA, infoB):

    posA = getpos(infoA)
    posB = getpos(infoB)
    if posB == '':
        posA, posB = getposend(infoA)
    return posA, posB
    

def get_homlen(info):
    homlen = 0
    l_info = info.split(';')
    for v in l_info:
        if v.startswith('HOMLEN='):
            homlen = v.replace('HOMLEN=','')
    return homlen


def get_homseq(info):
    homseq = ''
    l_info = info.split(';')
    for v in l_info:
        if v.startswith('HOMSEQ='):
            homseq = v.replace('HOMSEQ=','')
    return homseq


def sort_bp_bedpe(infile, outfile, file_type):
    
    hOUT = open(outfile, 'w')
    with open(infile, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            line = line.rstrip('\n')
            F = line.split('\t')

            chr1, start1, end1 = F[0], F[1], F[2]
            chr2, start2, end2 = F[3], F[4], F[5]
            dir1, dir2 = F[8], F[9]

            # The result of SvABA has error records.
            if "" in [chr1, start1, end1, dir1, chr2, start2, end2, dir2]: continue

            tmpChr1 = chr1 
            tmpChr2 = chr2
            if tmpChr1.startswith("chr"): tmpChr1 = tmpChr1.replace("chr","")
            if tmpChr2.startswith("chr"): tmpChr2 = tmpChr2.replace("chr","")
            if tmpChr1 == "X": tmpChr1 = "23"
            if tmpChr1 == "Y": tmpChr1 = "24"
            if tmpChr2 == "X": tmpChr2 = "23"
            if tmpChr2 == "Y": tmpChr2 = "24"
            
            val = line
            if (int(tmpChr1) > int(tmpChr2)) or (chr1 == chr2 and int(start1) > int(start2)):
                val = "\t".join([chr2, start2, end2, chr1, start1, end1, F[6], F[7], dir2, dir1])
                if file_type == "GenomonSV":
                    # GenomonSV
                    # 11.bp2_pos, 10.bp1_pos, 12.insert_seq, 13.sv_type, 15.bp2_gene, 14.bp1_gene, 17.bp2_exon, 16.bp1_exon, 28.bp2_overhung, 27,bp1_ovarhunb
                    val = val +"\t"+F[11]+"\t"+F[10]+"\t"+F[12]+"\t"+F[13]+"\t"+F[15]+"\t"+F[14]+"\t"+F[17]+"\t"+F[16]+"\t"+"\t".join(F[18:26]) +"\t"+F[28]+"\t"+F[27]
                elif file_type == "Manta" or file_type == "SvABA":
                    if F[10] == "BND":
                        # 10.sv_type, 11.Filter, 15.bp2.ID, 16.bp2.REF, 17.bp2_ALT, 12.bp1.ID, 13.bp1.REF, 14.bp1_ALT, 19.bp2_INFO, 28,bp1_INFO
                        val = val +"\t"+F[10]+"\t"+F[11]+"\t"+F[15]+"\t"+F[16]+"\t"+F[17]+"\t"+F[12]+"\t"+F[13]+"\t"+F[14]+"\t"+F[19]+"\t"+F[18]+"\t"+"\t".join(F[20:])
                    else:
                        val = val +"\t"+"\t".join(F[10:])
                else:
                    val = "file_type: " + file_type + " is not supported"
                    
            print(val, file=hOUT)


def sort_breakpoint(chr1,pos1,dir1,chr2,pos2,dir2):
    
    # The result of SvABA has error records.
    if "" in [chr1, pos1, dir1, chr2, pos2, dir2]:
        return chr1, pos1, dir1, chr2, pos2, dir2 
    
    tmpChr1 = chr1 
    tmpChr2 = chr2
    if tmpChr1.startswith("chr"): tmpChr1 = tmpChr1.replace("chr","")
    if tmpChr2.startswith("chr"): tmpChr2 = tmpChr2.replace("chr","")
    if tmpChr1 == "X": tmpChr1 = "23"
    if tmpChr1 == "Y": tmpChr1 = "24"
    if tmpChr2 == "X": tmpChr2 = "23"
    if tmpChr2 == "Y": tmpChr2 = "24"
    if (int(tmpChr1) > int(tmpChr2)) or (chr1 == chr2 and int(pos1) > int(pos2)):
        return chr2, pos2, dir2, chr1, pos1, dir1
    else:
        return chr1, pos1, dir1, chr2, pos2, dir2 


def get_hash_key(F, margin):
    
    chrA, startA, endA = F[0], str(int(F[1])+int(margin)), str(int(F[2])-int(margin))
    chrB, startB, endB = F[3], str(int(F[4])+int(margin)), str(int(F[5])-int(margin))
    strandA = F[8]
    strandB = F[9]
    # key 
    chrA, startA, strandA, chrB, startB, strandB = sort_breakpoint(chrA, startA, strandA, chrB, startB, strandB)
    key = chrA+'\t'+startA+'\t'+endA+'\t'+chrB+'\t'+startB+'\t'+endB+'\t'+strandA+'\t'+strandB
    return key


def get_hash_value(chrA, startA, strandA, chrB, startB, strandB, insert_seq):
    
    chrA, startA, strandA, chrB, startB, strandB = sort_breakpoint(chrA, startA, strandA, chrB, startB, strandB)
    val = chrA+','+startA+','+strandA+','+chrB+','+startB+','+strandB+','+insert_seq
    return val


def makeHash(anno_file, file_type_pair, margin, sv_type_idx, pair_idx, pair_info_idx, alt_idx):
    sv_comp = {}
    with open(anno_file, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')

            # get hash key
            key = get_hash_key(F, margin)
            
            # SV information (tool2)
            chrA_pair = F[pair_idx]
            startA_pair = str(int(F[pair_idx+1])+int(margin))
            chrB_pair = F[pair_idx+3]
            startB_pair = str(int(F[pair_idx+4])+int(margin))
            strandA_pair, strandB_pair, sv_type_pair = F[pair_idx+8], F[pair_idx+9], F[pair_idx+10]
            infoA_pair, infoB_pair = F[pair_info_idx], F[pair_info_idx+1]
            posA_pair, posB_pair = get_position(infoA_pair, infoB_pair)
            homlen = get_homlen(infoB_pair)
            homseq = get_homseq(infoB_pair)
            insert_seq = get_insert_seq(file_type_pair, F, pair_info_idx, alt_idx)

            val = ""
            # Manta
            if file_type_pair == "Manta":
                if sv_type_pair == "BND":
                    if strandA_pair == strandB_pair:
                        posB_pair = str(int(homlen) + int(posB_pair))
                else:
                    # Tandem Duplication
                    if strandA_pair == "-" and strandB_pair == "+": posA_pair = str(int(posA_pair)+1) 
                    # Deletion
                    if strandA_pair == "+" and strandB_pair == "-": posB_pair = str(int(posB_pair)+1)
                    
                val = get_hash_value(chrA_pair, posA_pair, strandA_pair, chrB_pair, posB_pair, strandB_pair, insert_seq)

            # SvABA
            elif file_type_pair == "SvABA":
                # The result of SvABA has error records.
                if posA_pair == "" or posB_pair == "": continue

                if infoB_pair == ".":
                    posB_pair = str(int(posB_pair)+1)
                else:
                    homlen = len(get_homseq(infoB_pair))
                    if strandB_pair == "+":
                        posB_pair = str(int(posB_pair) - int(homlen))
                    else:
                        posB_pair = str(int(posB_pair) + int(homlen))
                    
                val = get_hash_value(chrA_pair, posA_pair, strandA_pair, chrB_pair, posB_pair, strandB_pair, insert_seq)
                
            # GRIDSS
            elif file_type_pair == "GRIDSS":
                # posB_pair = str(int(homlen) + int(posB_pair))
    
                val = get_hash_value(chrA_pair, posA_pair, strandA_pair, chrB_pair, posB_pair, strandB_pair, insert_seq)
                
            # GenomonSV
            else:
                val = get_hash_value(chrA_pair, startA_pair, strandA_pair, chrB_pair, startB_pair, strandB_pair, insert_seq)
                
            sv_comp[key] = val
    return sv_comp


def annotate_other_sv(main_bedpe, outfile, sv_comp1, sv_comp2, sv_comp3, margin):
    
    hOUT = open(outfile, 'w')
    with open(main_bedpe) as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')
            # get hash key
            key = get_hash_key(F, margin)
            
            chrA, startA, endA = F[0], str(int(F[1])+int(margin)), str(int(F[2])-int(margin))
            chrB, startB, endB = F[3], str(int(F[4])+int(margin)), str(int(F[5])-int(margin))
            strandA = F[8]
            strandB = F[9]
            
            sv_info1 = sv_comp1[key] if key in sv_comp1 else "---"
            sv_info2 = sv_comp2[key] if key in sv_comp2 else "---"
            sv_info3 = sv_comp3[key] if key in sv_comp3 else "---"
            
            print( chrA+'\t'+startA+'\t'+endA+'\t'+chrB+'\t'+startB+'\t'+endB+'\t'+'\t'.join(F[6:]) +'\t'+ sv_info1 +'\t'+ sv_info2 +'\t'+ sv_info3, file=hOUT)
    hOUT.close()


# GenomonSV
def merge_bedpe1(annotate_bedpe, file_type, sv_hash):

    with open(annotate_bedpe, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')

            # get breakpint pair (chrA, posA, strandA, chrB, posB, strandB)
            chrA, chrB, posA, posB, strandA, strandB = F[0], F[3], F[10], F[11], F[8], F[9]
            insert_seq = get_insert_seq(file_type, F, 12, None)
            
            sv_type = get_sv_type(F[13])
            is_bedpe2 = F[29] # the bp pair of Manta
            is_bedpe3 =F[30] # the bp pair of SvABA
            is_bedpe4 =F[31] # the bp pair of GIRDSS

            chrA, posA, strandA, chrB, posB, strandB = sort_breakpoint(chrA, posA, strandA, chrB, posB, strandB)

            val = chrA+','+posA+','+strandA+','+chrB+','+posB+','+strandB+','+insert_seq # GenomonSV
            val = (val+"\t---") if is_bedpe2 == "---" else (val+"\t"+is_bedpe2) # Manta
            val = (val+"\t---") if is_bedpe3 == "---" else (val+"\t"+is_bedpe3) # SvABA
            val = (val+"\t---") if is_bedpe4 == "---" else (val+"\t"+is_bedpe4) # GRIDSS
            
            if insert_seq == '': insert_seq = '---'
            
            key = chrA+'\t'+posA+'\t'+strandA+'\t'+chrB+'\t'+posB+'\t'+strandB+'\t'+insert_seq+'\t'+sv_type

            
            sv_hash[key] = val
    return sv_hash

# Manta
def merge_bedpe2(annotate_bedpe, file_type, sv_hash, f_germ):

    with open(annotate_bedpe, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')

            # get breakpint pair (chrA, posA, strandA, chrB, posB, strandB)
            chrA, chrB, strandA, strandB = F[0], F[3], F[8], F[9]
            infoA, infoB = F[18], F[19]
            posA, posB = get_position(infoA, infoB)

            chrA, posA, strandA, chrB, posB, strandB = sort_breakpoint(chrA, posA, strandA, chrB, posB, strandB)

            sv_type = F[10]
            insert_seq = get_insert_seq(file_type, F, 18, None)
            if f_germ:
                is_bedpe1 = F[22] # the bp pair of GenomonSV
                is_bedpe3 =F[23] # the bp pair of SvABA
                is_bedpe4 =F[24] # the bp pair of Gridss
            else:
                is_bedpe1 = F[23] # the bp pair of GenomonSV
                is_bedpe3 =F[24] # the bp pair of SvABA
                is_bedpe4 =F[25] # the bp pair of GridSS
                
        
            # IF Genomon has already written this SV, continue
            if is_bedpe1 != "---": continue
        
            if sv_type == "BND":
                if strandA == strandB:
                    posB = str(int(get_homlen(infoB)) + int(posB))
            else:
                if strandA == "-" and strandB == "+": posA = str(int(posA)+1)
                if strandA == "+" and strandB == "-": posB = str(int(posB)+1)

            
            val = "---\t" # GenomonSV
            val = val + chrA+','+posA+','+strandA+','+chrB+','+posB+','+strandB+','+insert_seq # Manta
            val = (val+"\t---") if is_bedpe3 == "---" else (val+"\t"+is_bedpe3) # SvABA
            val = (val+"\t---") if is_bedpe4 == "---" else (val+"\t"+is_bedpe4) # GRIDSS
            
            if insert_seq == '': insert_seq = '---'
            key = chrA+'\t'+posA+'\t'+strandA+'\t'+chrB+'\t'+posB+'\t'+strandB+'\t'+insert_seq+'\t'+sv_type
            
            sv_hash[key] = val
    return sv_hash

# SvABA
def merge_bedpe3(annotate_bedpe, file_type, sv_hash, f_germ):

    with open(annotate_bedpe, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')
            
            # get breakpint pair (chrA, posA, strandA, chrB, posB, strandB)
            chrA, chrB, strandA, strandB = F[0], F[3], F[8], F[9]
            infoA, infoB = F[18], F[19]
            posA, posB = get_position(infoA, infoB)

            sv_type = F[10]
            insert_seq = get_insert_seq(file_type, F, 18, None)
            if f_germ:
                is_bedpe1 = F[22] # the bp pair of GenomonSV
                is_bedpe2 =F[23] # the bp pair of Manta
                is_bedpe4 =F[24] # the bp pair of GRIDSS
            else:
                is_bedpe1 = F[23] # the bp pair of GenomonSV
                is_bedpe2 =F[24] # the bp pair of Manta
                is_bedpe4 =F[25] # the bp pair of GRIDSS
            
            # the original result file of SvABA has some probrem.
            # skip the record of error.
            if posA == "" or posB == "": continue
            
            chrA, posA, strandA, chrB, posB, strandB = sort_breakpoint(chrA, posA, strandA, chrB, posB, strandB)
            
            # IF Genomon or Manta have already written this SV, continue
            if is_bedpe1 != "---" or is_bedpe2 != "---" : continue
        
            if infoB == ".":
                posB = str(int(posB)+1)
            else:
                homlen = len(get_homseq(infoB))
                if strandB == "+":
                    posB = str(int(posB) - int(homlen))
                else:
                    posB = str(int(posB) + int(homlen))

            val = "---\t" # GenomonSV
            val = val + "---\t" # Manta
            val = val + chrA+','+posA+','+strandA+','+chrB+','+posB+','+strandB+','+insert_seq # SvABA
            val = (val+"\t---") if is_bedpe4 == "---" else (val+"\t"+is_bedpe4) # GRIDSS

            if insert_seq == '': insert_seq = '---'
            key = chrA+'\t'+posA+'\t'+strandA+'\t'+chrB+'\t'+posB+'\t'+strandB+'\t'+insert_seq+'\t'+sv_type
            
            sv_hash[key] = val
    return sv_hash


# GRIDSS
def merge_bedpe4(annotate_bedpe, file_type, sv_hash, f_germ):

    with open(annotate_bedpe, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')

            # get breakpint pair (chrA, posA, strandA, chrB, posB, strandB)
            chrA, chrB, strandA, strandB = F[0], F[3], F[8], F[9]
            infoA, infoB = F[18], F[19]
            posA, posB = get_position(infoA, infoB)

            chrA, posA, strandA, chrB, posB, strandB = sort_breakpoint(chrA, posA, strandA, chrB, posB, strandB)

            sv_type = F[10]
            insert_seq = get_insert_seq(file_type, F, 18, 14)
            if f_germ:
                is_bedpe1 = F[22] # the bp pair of GenomonSV
                is_bedpe2 = F[23] # the bp pair of Manta
                is_bedpe3 =F[24] # the bp pair of SvABA
            else:
                is_bedpe1 = F[23] # the bp pair of GenomonSV
                is_bedpe2 = F[24] # the bp pair of Manta
                is_bedpe3 =F[25] # the bp pair of SvABA
                
            
            # IF Genomon, Manta or SvABA have already written this SV, continue
            if is_bedpe1 != "---" or is_bedpe2 != "---" or is_bedpe3 != "---" : continue
            
            val = "---\t" # GenomonSV
            val = val + "---\t" # Manta
            val = val + "---\t" # SvABA
            val = val + chrA+','+posA+','+strandA+','+chrB+','+posB+','+strandB+','+insert_seq # GRIDSS
            
            if insert_seq == '': insert_seq = '---'
            key = chrA+'\t'+posA+'\t'+strandA+'\t'+chrB+'\t'+posB+'\t'+strandB+'\t'+insert_seq+'\t'+sv_type
            
            sv_hash[key] = val
    return sv_hash



def simple_repeat_check(chr, start, end, simple_repeat_tb):

    sequence_array = []
    # check junction annotation for refGene for the first break point
    tabixErrorFlag = 0
    try:
        records = simple_repeat_tb.fetch(chr, int(start) - 1, int(end) + 1)
    except Exception as inst:
        print("%s: %s" % (type(inst), inst.args), file = sys.stderr)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if int(record[1]) <= int(start) + 1 and int(end) - 1 <= int(record[2]):
                sequence_array.append(record[15])
                
    return ";".join(sequence_array)


def get_simple_repeat(in_txt, output, simple_repeat_info):

    hOUT = open(output, 'w')
    with open(in_txt, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            
            chrA, chrB = F[0], F[3]
            posA, posB = F[1], F[4]
            simple_repeat_tb = pysam.TabixFile(simple_repeat_info)
            simpleA = simple_repeat_check(chrA, posA, posA, simple_repeat_tb)
            simpleB = simple_repeat_check(chrB, posB, posB, simple_repeat_tb)
            
            print(line + "\t"+ str(simpleA) + "\t"+ str(simpleB), file=hOUT)
    hOUT.close()

def print_result(in_txt, output):

    hOUT = open(output, 'w')
    print("Chr_1\tPos_1\tDir_1\tChr_2\tPos_2\tDir_2\tInserted_Seq\tVariant_Type\tis_Genomon\tis_Manta\tis_SvABA\tis_GRIDSS\tdist_to_exon\texon\tSimpleRepeat_1\tSimpleRepeat_2\tsv_size\tcount\tGene_1\tGene_2\tExon_1\tExon_2", file=hOUT)
    with open(in_txt, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            
            chrA, chrB = F[0], F[3]
            posA, posB = F[1], F[4]
            l_anot = [F[8], F[9], F[10], F[11]]
            size = int(posB) - int(posA) -1 if chrA == chrB and posA != '' and posB != '' else ''
            count = 4 - l_anot.count('---')
            
            print(line + "\t"+ str(size) + "\t"+ str(count) +"\t\t\t\t", file=hOUT)
    hOUT.close()


def merge_bp(in_sort_file, out_merged_file):
    
    hOUT = open(out_merged_file, 'w')
    with open(in_sort_file, 'r') as hIN:
        arch_bp = '' 
        arch_record = '' 
        for line in hIN:
            line = line.rstrip('\n')
            F = line.split('\t')

            bp = F[8]
            if bp == "---":
                print(line, file = hOUT)
                continue
            
            if arch_bp == '':
                arch_bp = bp
                arch_record = line
                
            elif check_bp(arch_bp, bp) == 0:
                # print("----------------------------------")
                # print(arch_bp)
                # print(bp)
                print(arch_record, file = hOUT)
                arch_bp = bp
                arch_record = line
                
            else: # if arch_bp and bp are same breakpoint position:
                print("----------------------------------")
                print(arch_bp)
                print(bp)
                arch_inseq = get_inseq_from_bp(arch_bp)
                inseq = get_inseq_from_bp(bp)
                if len(arch_inseq) > len(inseq):
                    arch_bp = bp
                    arch_record = line
                
        print(arch_record, file = hOUT)
    hOUT.close()    

def merge_SVs(args):

    print("START!")
    out_pref, ext = os.path.splitext(args.output)

    # sort break point TODO
    # sort_bp_bedpe(args.in_bedpe1, out_pref+".tmp1.sorted.bed", "GenomonSV")
    # sort_bp_bedpe(args.in_bedpe2, out_pref+".tmp2.sorted.bed", "Manta")
    # sort_bp_bedpe(args.in_bedpe3, out_pref+".tmp3.sorted.bed", "SvABA")

    hOUT = open(out_pref + ".tmp12.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe1, "-b", args.in_bedpe2], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp13.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe1, "-b", args.in_bedpe3], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp14.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe1, "-b", args.in_bedpe4], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp21.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe2, "-b", args.in_bedpe1], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp23.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe2, "-b", args.in_bedpe3], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp24.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe2, "-b", args.in_bedpe4], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp31.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe3, "-b", args.in_bedpe1], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp32.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe3, "-b", args.in_bedpe2], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp34.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe3, "-b", args.in_bedpe4], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp41.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe4, "-b", args.in_bedpe1], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp42.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe4, "-b", args.in_bedpe2], stdout = hOUT)
    hOUT.close()
    hOUT = open(out_pref + ".tmp43.bedpe", 'w')
    subprocess.check_call(["bedtools", "pairtopair", "-a", args.in_bedpe4, "-b", args.in_bedpe3], stdout = hOUT)
    hOUT.close()

    sv_comp1 = makeHash(out_pref + ".tmp12.bedpe", "Manta", args.margin, 13, 29, 47, None)
    sv_comp2 = makeHash(out_pref + ".tmp13.bedpe", "SvABA", args.margin, 13, 29, 47, None)
    sv_comp3 = makeHash(out_pref + ".tmp14.bedpe", "GRIDSS", args.margin, 13, 29, 47, 43)
    annotate_other_sv(args.in_bedpe1, out_pref + ".tmp1.annotate.bedpe", sv_comp1, sv_comp2, sv_comp3, args.margin)
    
    pair_idx = 22 if args.f_germ else 23
    info_gsv_idx = 34 if args.f_germ else 35
    info_bedpe_idx = 40 if args.f_germ else 41
    alt_idx = 37
    
    sv_comp1 = makeHash(out_pref + ".tmp21.bedpe", "GenomonSV", args.margin, 10, pair_idx, info_gsv_idx, None)
    sv_comp2 = makeHash(out_pref + ".tmp23.bedpe", "SvABA", args.margin, 10, pair_idx, info_bedpe_idx, None)
    sv_comp3 = makeHash(out_pref + ".tmp24.bedpe", "GRIDSS", args.margin, 10, pair_idx, info_bedpe_idx, alt_idx)
    annotate_other_sv(args.in_bedpe2, out_pref + ".tmp2.annotate.bedpe", sv_comp1, sv_comp2, sv_comp3, args.margin)

    sv_comp1 = makeHash(out_pref + ".tmp31.bedpe", "GenomonSV", args.margin, 10, pair_idx, info_gsv_idx, None)
    sv_comp2 = makeHash(out_pref + ".tmp32.bedpe", "Manta", args.margin, 10, pair_idx, info_bedpe_idx, None)
    sv_comp3 = makeHash(out_pref + ".tmp34.bedpe", "GRIDSS", args.margin, 10, pair_idx, info_bedpe_idx, alt_idx)
    annotate_other_sv(args.in_bedpe3, out_pref + ".tmp3.annotate.bedpe", sv_comp1, sv_comp2, sv_comp3, args.margin)

    sv_comp1 = makeHash(out_pref + ".tmp41.bedpe", "GenomonSV", args.margin, 10, pair_idx, info_gsv_idx, None)
    sv_comp2 = makeHash(out_pref + ".tmp42.bedpe", "Manta", args.margin, 10, pair_idx, info_bedpe_idx, None)
    sv_comp3 = makeHash(out_pref + ".tmp43.bedpe", "SvABA", args.margin, 10, pair_idx, info_bedpe_idx, None)
    annotate_other_sv(args.in_bedpe4, out_pref + ".tmp4.annotate.bedpe", sv_comp1, sv_comp2, sv_comp3, args.margin)

    sv_hash = {}
    sv_hash = merge_bedpe1(out_pref + ".tmp1.annotate.bedpe", "GenomonSV", sv_hash)
    sv_hash = merge_bedpe2(out_pref + ".tmp2.annotate.bedpe", "Manta", sv_hash, args.f_germ)
    sv_hash = merge_bedpe3(out_pref + ".tmp3.annotate.bedpe", "SvABA", sv_hash, args.f_germ)
    sv_hash = merge_bedpe4(out_pref + ".tmp4.annotate.bedpe", "GRIDSS", sv_hash, args.f_germ)
    
    
    with open(out_pref + ".tmp.annotate.txt", 'w') as hOUT:
        for key in sv_hash:
            print(key +'\t'+ sv_hash[key], file=hOUT)

    hOUT = open(out_pref + ".tmp.sorted.txt", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", out_pref + ".tmp.annotate.txt"],  stdout = hOUT)
    hOUT.close()
    
    merge_bp(out_pref + ".tmp.sorted.txt", out_pref + ".tmp.merged.txt")
            
    hOUT = open(out_pref + ".tmp.sorted2.txt", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", out_pref + ".tmp.merged.txt"],  stdout = hOUT)
    hOUT.close()

    cmd = []
    if args.genome_id == "hg19":
        if args.f_grc:
            cmd = ["sv_utils", "annotation", "--closest_exon", out_pref + ".tmp.sorted2.txt", out_pref + ".tmp.exon.txt"]
        else:
            cmd = ["sv_utils", "annotation", "--closest_exon", "--grc", out_pref + ".tmp.sorted2.txt", out_pref + ".tmp.exon.txt"]
    elif args.genome_id == "hg38":
        if args.f_grc:
            cmd = ["sv_utils", "annotation", "--genome_id", "hg38", "--closest_exon", out_pref + ".tmp.sorted2.txt", out_pref + ".tmp.exon.txt"]
        else:
            cmd = ["sv_utils", "annotation", "--genome_id", "hg38", "--closest_exon", "--grc", out_pref + ".tmp.sorted2.txt", out_pref + ".tmp.exon.txt"]
    subprocess.check_call(cmd)

    get_simple_repeat(out_pref + ".tmp.exon.txt", out_pref + ".tmp.simple.txt", args.simple_repeat_file)

    print_result(out_pref + ".tmp.simple.txt", out_pref + ".tmp.exon2.txt")
    
    if args.genome_id == "hg19":
        if args.f_grc:
            cmd = ["sv_utils", "annotation", "--re_gene_annotation", out_pref + ".tmp.exon2.txt", args.output]
        else:
            cmd = ["sv_utils", "annotation", "--re_gene_annotation", "--grc", out_pref + ".tmp.exon2.txt", args.output]
    elif args.genome_id == "hg38":
        if args.f_grc:
            cmd = ["sv_utils", "annotation", "--genome_id", "hg38", "--re_gene_annotation", out_pref + ".tmp.exon2.txt", args.output]
        else:
            cmd = ["sv_utils", "annotation", "--genome_id", "hg38", "--re_gene_annotation", "--grc", out_pref + ".tmp.exon2.txt", args.output]
    subprocess.check_call(cmd)



    os.remove(out_pref + ".tmp12.bedpe")
    os.remove(out_pref + ".tmp13.bedpe")
    os.remove(out_pref + ".tmp14.bedpe")
    os.remove(out_pref + ".tmp21.bedpe")
    os.remove(out_pref + ".tmp23.bedpe")
    os.remove(out_pref + ".tmp24.bedpe")
    os.remove(out_pref + ".tmp31.bedpe")
    os.remove(out_pref + ".tmp32.bedpe")
    os.remove(out_pref + ".tmp34.bedpe")
    os.remove(out_pref + ".tmp41.bedpe")
    os.remove(out_pref + ".tmp42.bedpe")
    os.remove(out_pref + ".tmp43.bedpe")
    os.remove(out_pref + ".tmp1.annotate.bedpe")
    os.remove(out_pref + ".tmp2.annotate.bedpe")
    os.remove(out_pref + ".tmp3.annotate.bedpe")
    os.remove(out_pref + ".tmp4.annotate.bedpe")
    os.remove(out_pref + ".tmp.annotate.txt")
    os.remove(out_pref + ".tmp.sorted.txt")
    os.remove(out_pref + ".tmp.merged.txt")
    # os.remove(out_pref + ".tmp.sorted2.txt")
    # os.remove(out_pref + ".tmp.exon.txt")
    
