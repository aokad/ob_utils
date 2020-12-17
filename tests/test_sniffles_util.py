#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import ob_utils
from ob_utils.utils import make_chrom_number_dict
from ob_utils.sniffles_util import snifflesSVtoBedpe
from ob_utils.sniffles_util import simplify_sniffles
from ob_utils.sniffles_util import filt_clustered_rearrangement2
from ob_utils.sniffles_util import repair_dup_strand

class TestSnifflesUtil_main(unittest.TestCase):

    def setUp(self):
        self.parser = ob_utils.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_vcf = cur_dir + "/data/test1_sniffles.txt"
        answer_file = cur_dir + "/data/test1_sniffles_answer.txt"
        output = tmp_dir + "/test1_sniffles.txt"
        snifflesSVtoBedpe(input_vcf, output, True, True, "PASS")
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_vcf = cur_dir + "/data/test2_sniffles.txt"
        answer_file = cur_dir + "/data/test2_sniffles_answer.txt"
        output = tmp_dir + "/test2_sniffles.txt"
        snifflesSVtoBedpe(input_vcf, output, True, True, "PASS")
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_sniffles.txt")
        input_bedpe = cur_dir + "/data/test3_sniffles.txt"
        answer_file = cur_dir + "/data/test3_sniffles_answer.txt"
        output = tmp_dir + "/test3_sniffles_out.txt"
        with open(output,"w") as hout:
            simplify_sniffles(input_bedpe, hout, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_sniffles.txt")
        input_bedpe = cur_dir + "/data/test4_sniffles.txt"
        answer_file = cur_dir + "/data/test4_sniffles_answer.txt"
        output = tmp_dir + "/test4_sniffles_out.txt"
        with open(output,"w") as hout:
            simplify_sniffles(input_bedpe, hout, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_sniffles.txt")
        input_bedpe = cur_dir + "/data/test5_sniffles.txt"
        answer_file = cur_dir + "/data/test5_sniffles_answer.txt"
        output = tmp_dir + "/test5_sniffles_out.txt"
        with open(output,"w") as hout:
            simplify_sniffles(input_bedpe, hout, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)    
        
    def test6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_sniffles.txt")
        input_tumor_bedpe = cur_dir + "/data/test6_sniffles_tumor.txt"
        input_simplify_control_bedpe = cur_dir + "/data/test6_sniffles_control.txt.gz"
        answer_file = cur_dir + "/data/test6_sniffles_answer.txt"
        output = tmp_dir + "test6_sniffles.txt"
        filt_clustered_rearrangement2(input_tumor_bedpe, output, input_simplify_control_bedpe, 0, 3, 1, 50, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_sniffles.txt")
        input_tumor_bedpe = cur_dir + "/data/test7_sniffles_tumor.txt"
        input_simplify_control_bedpe = cur_dir + "/data/test7_sniffles_control.txt.gz"
        answer_file = cur_dir + "/data/test7_sniffles_answer.txt"
        output = tmp_dir + "/test7_sniffles.txt"
        filt_clustered_rearrangement2(input_tumor_bedpe, output, input_simplify_control_bedpe, 1, 3, 1, 50, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test8(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_sniffles.txt")
        input_tumor_bedpe = cur_dir + "/data/test8_sniffles_tumor.txt"
        input_simplify_control_bedpe = cur_dir + "/data/test8_sniffles_control.txt.gz"
        answer_file = cur_dir + "/data/test8_sniffles_answer.txt"
        output = tmp_dir + "/test8_sniffles.txt"
        filt_clustered_rearrangement2(input_tumor_bedpe, output, input_simplify_control_bedpe, 50, 5, 1, 50, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)    
        
    def test9(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_sniffles.txt")
        input_tumor_bedpe = cur_dir + "/data/test9_sniffles_tumor.txt"
        input_simplify_control_bedpe = cur_dir + "/data/test9_sniffles_control.txt.gz"
        answer_file = cur_dir + "/data/test9_sniffles_answer.txt"
        output = tmp_dir + "/test9_sniffles.txt"
        filt_clustered_rearrangement2(input_tumor_bedpe, output, input_simplify_control_bedpe, 1, 3, 5, 50, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)   
        
    def test10(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_tumor_bedpe = cur_dir + "/data/test10_sniffles_tumor.txt"
        answer_file = cur_dir + "/data/test10_sniffles_answer.txt"
        output = tmp_dir + "/test10_sniffles.txt"
        repair_dup_strand(input_tumor_bedpe, output)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)   
