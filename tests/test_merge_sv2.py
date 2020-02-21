#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
from ob_utils import merge_sv2 as merge2


class TestMerge2SV(unittest.TestCase):

    def setUp(self):
        samtools = "/home/ubuntu/miniconda2/bin/samtools"

    def test1(self):
        bp_pair = "1,4669434,-,3,135016646,+,GG"
        chr1, pos1, dir1, chr2, pos2, dir2, inseq = merge2.get_position_from_bp(bp_pair)
        self.assertEqual("1",chr1)
        self.assertEqual("4669434",pos1)
        self.assertEqual("-",dir1)
        self.assertEqual("3",chr2)
        self.assertEqual("135016646",pos2)
        self.assertEqual("+",dir2)
        self.assertEqual("GG",inseq)
        
    def test2(self):
        bp_pair = "1,4669434,-,3,135016646,+,"
        inseq = merge2.get_inseq_from_bp(bp_pair)
        self.assertEqual("",inseq)
        
    def test3(self):
        bp_pair1 = "1,4669434,-,3,135016645,+,G"
        bp_pair2 = "1,4669434,-,3,135016646,+,"
        ret = merge2.check_bp(bp_pair1,bp_pair2)
        self.assertEqual(1,ret)
        
    def test4(self):
        bp_pair1 = "1,4669434,-,3,135016645,+,GG"
        bp_pair2 = "1,4669434,-,3,135016646,+,G"
        ret = merge2.check_bp(bp_pair1,bp_pair2)
        self.assertEqual(1,ret)
        
    def test5(self):
        bp_pair1 = "1,4669433,-,3,135016645,+,G"
        bp_pair2 = "1,4669434,-,3,135016646,+,G"
        ret = merge2.check_bp(bp_pair1,bp_pair2)
        self.assertEqual(1,ret)
        
    def test6(self):
        bp_pair1 = "1,4669433,-,2,135016645,+,G"
        bp_pair2 = "1,4669434,-,3,135016646,+,G"
        ret = merge2.check_bp(bp_pair1,bp_pair2)
        self.assertEqual(0,ret)
        
    def test7(self):
        infos = "POS=12345678"
        ret = merge2.get_info_val(infos,"POS")
        self.assertEqual('12345678',ret)

    def test8(self):
        infos = "POS=12345678;HOGEHOGE=fugafuga"
        ret = merge2.get_info_val(infos,"POS")
        self.assertEqual('12345678',ret)
        
    def test9(self):
        infos = "HOGEHOGE=fugafuga;POS=12345678"
        ret = merge2.get_info_val(infos,"POS")
        self.assertEqual('12345678',ret)
        
    def test10(self):
        F = ["AAAAAA","SVINSSEQ=CCCCCC;INSERTION=GGGGGG","[chr1:123456[ACAC","ACAC]chr1:123456]","]chr1:123456]TGTG","TGTG[chr1:123456["]
        ret = merge2.get_insert_seq("GenomonSV",F,0,None)
        self.assertEqual('AAAAAA',ret)
        ret = merge2.get_insert_seq("Manta",F,1,None)
        self.assertEqual('CCCCCC',ret)
        ret = merge2.get_insert_seq("SvABA",F,1,None)
        self.assertEqual('GGGGGG',ret)
        ret = merge2.get_insert_seq("GRIDSS",F,None,2)
        self.assertEqual('ACA',ret)
        ret = merge2.get_insert_seq("GRIDSS",F,None,3)
        self.assertEqual('CAC',ret)
        ret = merge2.get_insert_seq("GRIDSS",F,None,4)
        self.assertEqual('TGT',ret)
        ret = merge2.get_insert_seq("GRIDSS",F,None,5)
        self.assertEqual('GTG',ret)
        
    def test11(self):
        chr1, dir1, chr2, dir2 = "1","+","1","+"
        ret = merge2.get_sv_type_from_BND(chr1, dir1, chr2, dir2)
        self.assertEqual('INV',ret)
        
    def test12(self):
        chr1, dir1, chr2, dir2 = "1","+","1","-"
        ret = merge2.get_sv_type_from_BND(chr1, dir1, chr2, dir2)
        self.assertEqual('TRANS',ret)
        
    def test13(self):
        chr1, dir1, chr2, dir2 = "1","+","2","+"
        ret = merge2.get_sv_type_from_BND(chr1, dir1, chr2, dir2)
        self.assertEqual('TRANS',ret)
        
    def test14(self):
        info = "POS=1111111;END=222222"
        pos, end = merge2.getposend(info)
        self.assertEqual('1111111',pos)
        self.assertEqual('222222',end)
        
    def test15(self):
        info = ""
        pos, end = merge2.getposend(info)
        self.assertEqual('',pos)
        self.assertEqual('',end)
        
    def test16(self):
        info = "POS=1111111;END=222222;HOMLEN=333333"
        homlen = merge2.get_homlen(info)
        self.assertEqual('333333',homlen)
        
    def test17(self):
        info = "POS=1111111;END=222222;HOMLEN=333333;HOMSEQ=444444"
        homseq = merge2.get_homseq(info)
        self.assertEqual('444444',homseq)

    def test18(self):
        bp_pair = "Y,4669434,-,X,135016646,+,GG"
        chr1, pos1, dir1, chr2, pos2, dir2, inseq = merge2.get_position_from_bp(bp_pair)
        chr1, pos1, dir1, chr2, pos2, dir2 = merge2.sort_breakpoint(chr1, pos1, dir1, chr2, pos2, dir2)
        self.assertEqual("X",chr1)
        self.assertEqual("135016646",pos1)
        self.assertEqual("+",dir1)
        self.assertEqual("Y",chr2)
        self.assertEqual("4669434",pos2)
        self.assertEqual("-",dir2)
        self.assertEqual("GG",inseq)

    def test19(self):
        bp_pair = "1,4669434,-,X,135016646,+,GG"
        chr1, pos1, dir1, chr2, pos2, dir2, inseq = merge2.get_position_from_bp(bp_pair)
        chr1, pos1, dir1, chr2, pos2, dir2 = merge2.sort_breakpoint(chr1, pos1, dir1, chr2, pos2, dir2)
        self.assertEqual("1",chr1)
        self.assertEqual("4669434",pos1)
        self.assertEqual("-",dir1)
        self.assertEqual("X",chr2)
        self.assertEqual("135016646",pos2)
        self.assertEqual("+",dir2)
        self.assertEqual("GG",inseq)


    '''
    def test20(self):
        bp_pair = "1,4669434,-,X,135016646,+,GG"

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        anno_file = cur_dir + "/../data/COLO829.merged.somaticSV.tmp12.bedpe"
        file_type_pair = "Manta"
        margin = 10
        sv_type_idx = 13
        pair_idx = 29
        pair_info_idx = 47
        sv_comp = merge2.makeHash(anno_file, file_type_pair, margin, sv_type_idx, pair_idx, pair_info_idx)
        ret = sv_comp(1	207981223	207981243	1	208014812	208014832)
        self.assertEqual("1,207981229,+,1,208014817,-",ret)
        
    '''    