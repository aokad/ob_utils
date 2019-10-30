#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import ob_utils


class TestCompSV(unittest.TestCase):

    def setUp(self):
        self.parser = ob_utils.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt1.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt1.txt"
        output_file = tmp_dir+"/5929_small_tumor_result1.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result1_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt2.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt2.txt"
        output_file = tmp_dir+"/5929_small_tumor_result2.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result2_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt3.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt3.txt"
        output_file = tmp_dir+"/5929_small_tumor_result3.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result3_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt4.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt4.txt"
        output_file = tmp_dir+"/5929_small_tumor_result4.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result4_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt5.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt5.txt"
        output_file = tmp_dir+"/5929_small_tumor_result5.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result5_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt6.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt6.txt"
        output_file = tmp_dir+"/5929_small_tumor_result6.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result6_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt7.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt7.txt"
        output_file = tmp_dir+"/5929_small_tumor_result7.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result7_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test8(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt8.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt8.txt"
        output_file = tmp_dir+"/5929_small_tumor_result8.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result8_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test9(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt9.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt9.txt"
        output_file = tmp_dir+"/5929_small_tumor_result9.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result9_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test10(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak = cur_dir + "/../data/5929_small_tumor_filt.txt"
        in_genomonsv = cur_dir + "/../data/5929_tumor.genomonSV.result.filt.txt"
        output_file = tmp_dir+"/5929_small_tumor_result.txt"
        args = self.parser.parse_args(["comp", "--in_onebreak", in_onebreak, "--in_genomonsv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_result_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        