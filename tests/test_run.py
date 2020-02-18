#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import ob_utils


class TestMantaSVtoBedpe_main(unittest.TestCase):

    def setUp(self):
        self.parser = ob_utils.parser.create_parser()

    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_manta_sv = cur_dir + "/../data/COLO829.manta.somtaicSV_test1.vcf"
        output_file = tmp_dir+"/COLO829.manta.somtaicSV_test1_out.bedpe"
        args = self.parser.parse_args(["manta_sv", "--in_manta_sv", in_manta_sv, "--output", output_file])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.manta.somtaicSV_test1_answer.bedpe"
        # self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_manta_sv = cur_dir + "/../data/COLO829.manta.somtaicSV_test1.vcf"
        output_file = tmp_dir+"/COLO829.manta.somtaicSV_test2_out.bedpe"
        args = self.parser.parse_args(["manta_sv", "--in_manta_sv", in_manta_sv, "--output", output_file, "--margin", "1"])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.manta.somtaicSV_test2_answer.bedpe"
        # self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_manta_sv = cur_dir + "/../data/COLO829.manta.somtaicSV_test2.vcf"
        output_file = tmp_dir+"/COLO829.manta.somtaicSV_test3_out.bedpe"
        args = self.parser.parse_args(["manta_sv", "--in_manta_sv", in_manta_sv, "--output", output_file, "--margin", "1"])
        args.func(args)
        answer_file = cur_dir + "/../data/COLO829.manta.somtaicSV_test3_answer.bedpe"
        # self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_manta_sv = cur_dir + "/../data/COLO829.manta.somtaicSV_test4.vcf"
        output_file = tmp_dir+"/COLO829.manta.somtaicSV_test4_out.bedpe"
        args = self.parser.parse_args(["manta_sv", "--in_manta_sv", in_manta_sv, "--output", output_file, "--margin", "1"])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.manta.somtaicSV_test4_answer.bedpe"
        # self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)


    def test5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_manta_sv = cur_dir + "/../data/COLO829.manta.somtaicSV_test5.vcf"
        output_file = tmp_dir+"/COLO829.manta.somtaicSV_test5_out.bedpe"
        args = self.parser.parse_args(["manta_sv", "--in_manta_sv", in_manta_sv, "--output", output_file, "--margin", "1"])
        args.func(args)

        answer_file = cur_dir + "/../data/COLO829.manta.somtaicSV_test5_answer.bedpe"
        # self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)
        
        