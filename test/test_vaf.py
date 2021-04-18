import os
import pysam
from nose.tools import *
from case_control_filter.vaf import get_vaf_method, _get_ad_vaf
from case_control_filter.vaf import _get_platypus_vaf, _get_svaba_vaf
from case_control_filter.vaf import _get_strelka_vaf, _get_freebayes_vaf

dir_path = os.path.dirname(os.path.realpath(__file__))
ad_vcf = os.path.join(dir_path, 'test_data', 'ad_test.vcf')
nv_vcf = os.path.join(dir_path, 'test_data', 'nv_test.vcf')
fb_vcf = os.path.join(dir_path, 'test_data', 'fb_test.vcf')
svaba_vcf = os.path.join(dir_path, 'test_data', 'svaba_test.vcf')
strelka_vcf = os.path.join(dir_path, 'test_data', 'strelka_test.vcf')


def test_get_ad_vaf():
    with pysam.VariantFile(ad_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_ad_vaf)


def test_get_nv_vaf():
    with pysam.VariantFile(nv_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_platypus_vaf)


def test_get_freebayes_vaf():
    with pysam.VariantFile(fb_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_freebayes_vaf)


def test_get_svaba_vaf():
    with pysam.VariantFile(svaba_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_svaba_vaf)


def test_get_strelka_vaf():
    with pysam.VariantFile(strelka_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_strelka_vaf)
