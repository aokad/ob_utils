#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: ken0-1n
"""

import sys
import argparse
from .version import __version__
from .comp_sv import comp_main
from .merge_sv import merge_main
from .genomonsv_util import genomonSVtoBedpe_main
from .manta_util import mantaSVtoBedpe_main
from .gridss_util import gridssSVtoBedpe_main
from .svaba_util import svabaSVtoBedpe_main
from .sniffles_util import snifflesSVtoBedpe_main
from .delly_util import dellySVtoBedpe_main
from .cutesv_util import cutesvSVtoBedpe_main
from .camphor_util import camphorSVtoBedpe_main
from .svim_util import svimSVtoBedpe_main
from .savana_util import savanaSVtoBedpe_main
from .savana103_util import savanaSVtoBedpe_main
from .common_bedpe_util import bedpetoBedpe_main
from .merge_sv2 import merge_SVs
from .liftover_trafic import liftover_trafic_main

def create_parser():
    prog = "ob_utils"
    parser = argparse.ArgumentParser(prog = prog)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()
    
    def _create_comp_parser(subparsers):
        
        comp_parser = subparsers.add_parser("comp", help = "onebreak results and GenomonSV results")
        comp_parser.add_argument("--in_onebreak", help = "the result of onebreak", type = str, required=True)
        comp_parser.add_argument("--in_genomonsv", help = "the result of GenomonSV", type = str, required=True)
        comp_parser.add_argument("--output", help = "the output file", type = str, required=True)
        comp_parser.add_argument("--margin", help = "the margin for comparing SVs and SVs", type = int, default = 10)
        return comp_parser

    
    def _create_merge_parser(subparsers):
        
        merge_parser = subparsers.add_parser("merge", help = "merge onebreak filt results")
        merge_parser.add_argument("--in_onebreak_filt1", help = "the result of onebreak filt 1", type = str, required=True)
        merge_parser.add_argument("--in_onebreak_filt2", help = "the result of onebreak filt 2", type = str, required=True)
        merge_parser.add_argument("--output", help = "tihe output file", type = str, required=True)
        return merge_parser
        

    def _create_genomonsv_util_parser(subparsers):
        
        genomonsv_parser = subparsers.add_parser("genomon_sv", help = "convert the genomon sv format to BEDPE file")
        genomonsv_parser.add_argument("--in_genomon_sv", help = "the bedpe format file", type = str, required=True)
        genomonsv_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        genomonsv_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 10)
        genomonsv_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        genomonsv_parser.add_argument("--v2", help = 'chromosome of sv file. True=new format|False=old format', action = 'store_true', default = False )        
        return genomonsv_parser
    
    
    def _create_manta_util_parser(subparsers):
        
        manta_parser = subparsers.add_parser("manta_sv", help = "convert the manta sv format to BEDPE file")
        manta_parser.add_argument("--in_manta_sv", help = "the bedpe format file", type = str, required=True)
        manta_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        manta_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 10)
        manta_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        manta_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        manta_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        manta_parser.add_argument("--is_cnv", default = False, action = 'store_true', help = "analysis CNV")
        return manta_parser
    
    
    def _create_gridss_util_parser(subparsers):
        
        gridss_parser = subparsers.add_parser("gridss_sv", help = "convert the manta sv format to BEDPE file")
        gridss_parser.add_argument("--in_gridss_sv", help = "the bedpe format file", type = str, required=True)
        gridss_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        gridss_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 10)
        gridss_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        gridss_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        gridss_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")
        gridss_parser.add_argument("--singleend_sv", default = False, action = 'store_true', help = "analysis single end SV")
        gridss_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        return gridss_parser
        
    
    def _create_svaba_util_parser(subparsers):
        
        svaba_parser = subparsers.add_parser("svaba_sv", help = "convert the svaba sv format to BEDPE file")
        svaba_parser.add_argument("--in_svaba_sv", help = "the bedpe format file", type = str, required=True)
        svaba_parser.add_argument("--in_svaba_indel", help = "the bedpe format file", type = str, required=True)
        svaba_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        svaba_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 10)
        svaba_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )
        svaba_parser.add_argument("--f_germ", help = 'sample is not Tumor/Normal pair.', action = 'store_true', default = False )
        svaba_parser.add_argument("--normal_max_variant", help = "max variant reads", type = int, default = 1)
        svaba_parser.add_argument("--tumor_min_variant", help = "min variant reads", type = int, default = 2)
        svaba_parser.add_argument("--normal_min_depth", help = "tumor depth", type = int, default = 10)
        svaba_parser.add_argument("--tumor_min_depth", help = "normal_depth", type = int, default = 10)
        svaba_parser.add_argument("--min_del_size", help = "min deleltion size", type = int, default = 20)
        svaba_parser.add_argument("--min_ins_size", help = "min insertion size", type = int, default = 12)
        return svaba_parser


    def _create_sniffles_util_parser(subparsers):
        
        sniffles_parser = subparsers.add_parser("sniffles_sv", help = "convert the sniffles sv format to BEDPE file")
        sniffles_parser.add_argument("--in_sniffles_tumor_sv", help = "the vcf format file", type = str, required=True)
        sniffles_parser.add_argument("--in_sniffles_control_sv", help = "the vcf format file", type = str, required=True)
        sniffles_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        sniffles_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 50)
        sniffles_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        sniffles_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        sniffles_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        sniffles_parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
        sniffles_parser.add_argument("--max_control_support_read", help = "maximum control support reads", type = int, default = 1)
        sniffles_parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
        sniffles_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")
        sniffles_parser.add_argument("--sniffles2", default = False, action = 'store_true', help = "used sniffles2")

        return sniffles_parser

    def _create_delly_util_parser(subparsers):
        
        delly_parser = subparsers.add_parser("delly_sv", help = "convert the delly sv format to BEDPE file")
        delly_parser.add_argument("--in_delly_tumor_sv", help = "the vcf format file", type = str, required=True)
        delly_parser.add_argument("--in_delly_control_sv", help = "the vcf format file", type = str, required=True)
        delly_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        delly_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 50)
        delly_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        delly_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        delly_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        delly_parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
        delly_parser.add_argument("--max_control_support_read", help = "maximum control support reads", type = int, default = 1)
        delly_parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
        delly_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

        return delly_parser

    def _create_cutesv_util_parser(subparsers):
        
        cutesv_parser = subparsers.add_parser("cutesv_sv", help = "convert the cuteSV sv format to BEDPE file")
        cutesv_parser.add_argument("--in_cutesv_tumor_sv", help = "the vcf format file", type = str, required=True)
        cutesv_parser.add_argument("--in_cutesv_control_sv", help = "the vcf format file", type = str, required=True)
        cutesv_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        cutesv_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 50)
        cutesv_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        cutesv_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        cutesv_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        cutesv_parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
        cutesv_parser.add_argument("--max_control_support_read", help = "maximum control support reads", type = int, default = 1)
        cutesv_parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
        cutesv_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

        return cutesv_parser

    def _create_camphor_util_parser(subparsers):
        
        camphor_parser = subparsers.add_parser("camphor_sv", help = "convert the CAMPHORsomatic sv format to BEDPE file")
        camphor_parser.add_argument("--in_camphor_tumor_sv", help = "the vcf format file", type = str, required=True)
        camphor_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        camphor_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        camphor_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        camphor_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        camphor_parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
        camphor_parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
        camphor_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

        return camphor_parser

    def _create_svim_util_parser(subparsers):
        
        svim_parser = subparsers.add_parser("svim_sv", help = "convert the svim sv format to BEDPE file")
        svim_parser.add_argument("--in_svim_tumor_sv", help = "the vcf format file", type = str, required=True)
        svim_parser.add_argument("--in_svim_control_sv", help = "the vcf format file", type = str, required=True)
        svim_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        svim_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 50)
        svim_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        svim_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        svim_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        svim_parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
        svim_parser.add_argument("--max_control_support_read", help = "maximum control support reads", type = int, default = 1)
        svim_parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
        svim_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

        return svim_parser

    def _create_savana_util_parser(subparsers):
        
        savana_parser = subparsers.add_parser("savana_sv", help = "convert the savana sv format to BEDPE file")
        savana_parser.add_argument("--in_savana_tumor_sv", help = "the vcf format file", type = str, required=True)
        savana_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        savana_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        savana_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        savana_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        savana_parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
        savana_parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
        savana_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

        return savana_parser

    def _create_savana103_util_parser(subparsers):
        
        savana_parser = subparsers.add_parser("savana1.0.3_sv", help = "convert the savana sv format to BEDPE file")
        savana_parser.add_argument("--in_savana_tumor_sv", help = "the vcf format file", type = str, required=True)
        savana_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        savana_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False )        
        savana_parser.add_argument("--bcf_filter_option", help = "filter options for bcftools view", type = str, default = "PASS")
        savana_parser.add_argument("--filter_scaffold_option", default = False, action = 'store_true', help = "if True, output only chr1-22 and XY.")
        savana_parser.add_argument("--min_tumor_support_read", help = "minimum tumor support reads", type = int, default = 3)
        savana_parser.add_argument("--min_sv_length", help = "minimum sv length", type = int, default = 1)
        savana_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

        return savana_parser

    def _create_common_bedpe_util_parser(subparsers):
        
        common_bedpe_parser = subparsers.add_parser("common_bedpe", help = "convert the common bedpe to BEDPE file for comparison")
        common_bedpe_parser.add_argument("--in_bedpe", help = "the bedpe format file", type = str, required=True)
        common_bedpe_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        common_bedpe_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 10)
        common_bedpe_parser.add_argument("--min_sv_size", help = "minimum sv size", type = int, default = 50)
        return common_bedpe_parser
    
    
    def _merge_sv_parser(subparsers):
        
        merge_svs_parser = subparsers.add_parser("merge_sv", help = "merge bedpe files")
        merge_svs_parser.add_argument("--in_genomonsv", help = "the bedpe format file (GenomonSV)", type = str, required=True)
        merge_svs_parser.add_argument("--in_manta", help = "the bedpe format file (Manta)", type = str, required=True)
        merge_svs_parser.add_argument("--in_svaba", help = "the bedpe format file (SvABA)", type = str, required=True)
        merge_svs_parser.add_argument("--in_gridss", help = "the bedpe format file (GRIDSS)", type = str, required=True)
        merge_svs_parser.add_argument("--output", help = "the output bedpe format file", type = str, required=True)
        merge_svs_parser.add_argument("--margin", help = "the margin for Bedpe", type = int, default = 10)
        merge_svs_parser.add_argument("--f_grc", help = 'chromosome of sv file. True=chr1|False=1', action = 'store_true', default = False ) 
        merge_svs_parser.add_argument("--f_germ", help = 'sample is not Tumor/Normal pair.', action = 'store_true', default = False )
        merge_svs_parser.add_argument("--simple_repeat_file", help = ' variants with overlapping simple repeat annotatio',  default = None, type = str , required=True)
        merge_svs_parser.add_argument("--reference", help = 'the reference genome',  default = None, type = str , required=True)
        merge_svs_parser.add_argument("--genome_id", help = 'the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: hg38).', choices = ["hg19", "hg38"], default = "hg38" )
        return merge_svs_parser
    
    
    def _liftover_trafic_parser(subparsers):
        
        lo_trafic_parser = subparsers.add_parser("liftover_trafic", help = "merge bedpe files")
        lo_trafic_parser.add_argument("--in_vcf", help = "the result of trafic (hg19)", type = str, required=True)
        lo_trafic_parser.add_argument("--out_vcf", help = 'the result of trafic (liftover hg38)', type = str, required=True) 
        lo_trafic_parser.add_argument("--map_chain", help = 'the map chain file', type = str, required=True, default = "hg19ToHg38.over.chain")
        lo_trafic_parser.add_argument("--debug", help = 'if True, not remove temp files', action = 'store_true', default = False )
        return lo_trafic_parser
    

    comp_parser = _create_comp_parser(subparsers)
    comp_parser.set_defaults(func = comp_main)
    merge_parser = _create_merge_parser(subparsers)
    merge_parser.set_defaults(func = merge_main)
    genomonsv_parser = _create_genomonsv_util_parser(subparsers)
    genomonsv_parser.set_defaults(func = genomonSVtoBedpe_main)
    manta_parser = _create_manta_util_parser(subparsers)
    manta_parser.set_defaults(func = mantaSVtoBedpe_main)
    gridss_parser = _create_gridss_util_parser(subparsers)
    gridss_parser.set_defaults(func = gridssSVtoBedpe_main)
    svaba_parser = _create_svaba_util_parser(subparsers)
    svaba_parser.set_defaults(func = svabaSVtoBedpe_main)
    sniffles_parser = _create_sniffles_util_parser(subparsers)
    sniffles_parser.set_defaults(func = snifflesSVtoBedpe_main)
    delly_parser = _create_delly_util_parser(subparsers)
    delly_parser.set_defaults(func = dellySVtoBedpe_main)
    cutesv_parser = _create_cutesv_util_parser(subparsers)
    cutesv_parser.set_defaults(func = cutesvSVtoBedpe_main)
    camphor_parser = _create_camphor_util_parser(subparsers)
    camphor_parser.set_defaults(func = camphorSVtoBedpe_main)
    common_bedpe_parser = _create_common_bedpe_util_parser(subparsers)
    common_bedpe_parser.set_defaults(func = bedpetoBedpe_main)
    svim_parser = _create_svim_util_parser(subparsers)
    svim_parser.set_defaults(func = svimSVtoBedpe_main)
    savana_parser = _create_savana_util_parser(subparsers)
    savana_parser.set_defaults(func = savanaSVtoBedpe_main)
    savana103_parser = _create_savana103_util_parser(subparsers)
    savana103_parser.set_defaults(func = savanaSVtoBedpe_main)
    merge_svs_parser = _merge_sv_parser(subparsers)
    merge_svs_parser.set_defaults(func = merge_SVs)
    lo_trafic_parser = _liftover_trafic_parser(subparsers)
    lo_trafic_parser.set_defaults(func = liftover_trafic_main)
    return parser
