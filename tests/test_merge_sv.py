#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import ob_utils


class TestMergeSV(unittest.TestCase):

    def setUp(self):
        self.parser = ob_utils.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak1 = cur_dir + "/../data/5929_small_tumor_filt1-1.txt"
        in_onebreak2 = cur_dir + "/../data/5929_small_tumor_filt1-2.txt"
        output_file = tmp_dir+"/5929_small_tumor_filt1-1_result.txt"
        args = self.parser.parse_args(["merge", "--in_onebreak_filt1", in_onebreak1, "--in_onebreak_filt2", in_onebreak2, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_1-1_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak1 = cur_dir + "/../data/5929_small_tumor_filt2-1.txt"
        in_onebreak2 = cur_dir + "/../data/5929_small_tumor_filt2-2.txt"
        output_file = tmp_dir+"/5929_small_tumor_filt2-1_result.txt"
        args = self.parser.parse_args(["merge", "--in_onebreak_filt1", in_onebreak1, "--in_onebreak_filt2", in_onebreak2, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_2-1_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak1 = cur_dir + "/../data/5929_small_tumor_filt3-1.txt"
        in_onebreak2 = cur_dir + "/../data/5929_small_tumor_filt3-2.txt"
        output_file = tmp_dir+"/5929_small_tumor_filt3-1_result.txt"
        args = self.parser.parse_args(["merge", "--in_onebreak_filt1", in_onebreak1, "--in_onebreak_filt2", in_onebreak2, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_3-1_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_onebreak1 = cur_dir + "/../data/5929_small_tumor_filt4-1.txt"
        in_onebreak2 = cur_dir + "/../data/5929_small_tumor_filt4-2.txt"
        output_file = tmp_dir+"/5929_small_tumor_filt4-1_result.txt"
        args = self.parser.parse_args(["merge", "--in_onebreak_filt1", in_onebreak1, "--in_onebreak_filt2", in_onebreak2, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_tumor_4-1_answer.txt"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

