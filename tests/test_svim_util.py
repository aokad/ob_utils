#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import ob_utils
from ob_utils.utils import make_chrom_number_dict
from ob_utils.svim_util import svimSVtoBedpe
from ob_utils.svim_util import simplify_svim
from ob_utils.svim_util import filt_clustered_rearrangement2

class TestSvimUtil_main(unittest.TestCase):

    def setUp(self):
        self.parser = ob_utils.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_vcf = cur_dir + "/data/test1_svim.txt"
        answer_file = cur_dir + "/data/test1_svim_answer.txt"
        output = tmp_dir + "/test1_svim.txt"
        svimSVtoBedpe(input_vcf, output, True, True, "PASS")
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        input_vcf = cur_dir + "/data/test2_svim.txt"
        answer_file = cur_dir + "/data/test2_svim_answer.txt"
        output = tmp_dir + "/test2_svim.txt"
        svimSVtoBedpe(input_vcf, output, True, True, "PASS")
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_svim.txt")
        input_bedpe = cur_dir + "/data/test3_svim.txt"
        answer_file = cur_dir + "/data/test3_svim_answer.txt"
        output = tmp_dir + "/test3_svim.txt"
        with open(output,"w") as hout:
            simplify_svim(input_bedpe, 1, hout, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_svim.txt")
        input_bedpe = cur_dir + "/data/test4_svim.txt"
        answer_file = cur_dir + "/data/test4_svim_answer.txt"
        output = tmp_dir + "/test4_svim.txt"
        with open(output,"w") as hout:
            simplify_svim(input_bedpe, 4, hout, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_svim.txt")
        input_bedpe = cur_dir + "/data/test5_svim.txt"
        answer_file = cur_dir + "/data/test5_svim_answer.txt"
        output = tmp_dir + "/test5_svim.txt"
        with open(output,"w") as hout:
            simplify_svim(input_bedpe, 1, hout, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_svim.txt")
        input_tumor_bedpe = cur_dir + "/data/test6_svim_tumor.txt"
        input_simplify_control_bedpe = cur_dir + "/data/test6_svim_control.txt.gz"
        answer_file = cur_dir + "/data/test6_svim_answer.txt"
        output = tmp_dir + "/test6_svim.txt"
        filt_clustered_rearrangement2(input_tumor_bedpe, output, input_simplify_control_bedpe, 0, 3, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_svim.txt")
        input_tumor_bedpe = cur_dir + "/data/test7_svim_tumor.txt"
        input_simplify_control_bedpe = cur_dir + "/data/test7_svim_control.txt.gz"
        answer_file = cur_dir + "/data/test7_svim_answer.txt"
        output = tmp_dir + "/test7_svim.txt"
        filt_clustered_rearrangement2(input_tumor_bedpe, output, input_simplify_control_bedpe, 0, 3, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test8(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        h_chrom_number = make_chrom_number_dict(cur_dir + "/data/test1_svim.txt")
        input_tumor_bedpe = cur_dir + "/data/test8_svim_tumor.txt"
        input_simplify_control_bedpe = cur_dir + "/data/test8_svim_control.txt.gz"
        answer_file = cur_dir + "/data/test8_svim_answer.txt"
        output = tmp_dir + "/test8_svim.txt"
        filt_clustered_rearrangement2(input_tumor_bedpe, output, input_simplify_control_bedpe, 0, 4, h_chrom_number)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
