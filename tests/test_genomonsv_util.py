#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import ob_utils


class TestGenomonSVtoBedpe_main(unittest.TestCase):

    def setUp(self):
        self.parser = ob_utils.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_genomonsv = cur_dir + "/../data/COLO829.genomonSV.result.filt.test1.txt"
        output_file = tmp_dir+"/COLO829.genomonSV.result.filt.test1_out.bedpe"
        args = self.parser.parse_args(["genomon_sv", "--in_genomon_sv", in_genomonsv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.genomonSV.result.filt.test1_answer.bedpe"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_genomonsv = cur_dir + "/../data/COLO829.genomonSV.result.filt.test1.txt"
        output_file = tmp_dir+"/COLO829.genomonSV.result.filt.test2_out.bedpe"
        args = self.parser.parse_args(["genomon_sv", "--in_genomon_sv", in_genomonsv, "--output", output_file, "--margin", "1"])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.genomonSV.result.filt.test2_answer.bedpe"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_genomonsv = cur_dir + "/../data/COLO829.genomonSV.result.filt.test2.txt"
        output_file = tmp_dir+"/COLO829.genomonSV.result.filt.test3_out.bedpe"
        args = self.parser.parse_args(["genomon_sv", "--in_genomon_sv", in_genomonsv, "--output", output_file, "--margin", "1"])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.genomonSV.result.filt.test3_answer.bedpe"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_genomonsv = cur_dir + "/../data/COLO829.genomonSV.result.filt.test4.txt"
        output_file = tmp_dir+"/COLO829.genomonSV.result.filt.test4_out.bedpe"
        args = self.parser.parse_args(["genomon_sv", "--in_genomon_sv", in_genomonsv, "--output", output_file, "--margin", "1"])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.genomonSV.result.filt.test4_answer.bedpe"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
