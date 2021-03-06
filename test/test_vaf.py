import os
import pysam
from nose.tools import *
from .utils import get_variants
from cafex.vaf import get_vaf_method, _get_ad_vaf
from cafex.vaf import _get_platypus_vaf, _get_svaba_vaf
from cafex.vaf import _get_strelka_vaf, _get_freebayes_vaf

dir_path = os.path.dirname(os.path.realpath(__file__))
ad_vcf = os.path.join(dir_path, 'test_data', 'ad_test.vcf')
nv_vcf = os.path.join(dir_path, 'test_data', 'nv_test.vcf')
fb_vcf = os.path.join(dir_path, 'test_data', 'fb_test.vcf')
svaba_vcf = os.path.join(dir_path, 'test_data', 'svaba_test.vcf')
strelka_vcf = os.path.join(dir_path, 'test_data', 'strelka_test.vcf')
no_vaf_vcf = os.path.join(dir_path, 'test_data', 'no_vaf_test.header.vcf')


def test_no_vaf_field():
    with pysam.VariantFile(no_vaf_vcf) as vcf:
        assert_raises(ValueError, get_vaf_method, vcf)
    

def test_get_ad_vaf_method():
    with pysam.VariantFile(ad_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_ad_vaf)


def test_get_nv_vaf_method():
    with pysam.VariantFile(nv_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_platypus_vaf)


def test_get_freebayes_vaf_method():
    with pysam.VariantFile(fb_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_freebayes_vaf)


def test_get_svaba_vaf_method():
    with pysam.VariantFile(svaba_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_svaba_vaf)


def test_get_strelka_vaf_method():
    with pysam.VariantFile(strelka_vcf) as vcf:
        result = get_vaf_method(vcf)
        assert_equal(result, _get_strelka_vaf)


def test_ad_vaf_values():
    ''' VAF from AD '''
    expected = {
        'Case1': [[0.4], [0.4666666666666667], [1.0], [0.5], [1.0], [0.0, 0.5],
                  [0.21875], [1.0], [0.0], [0.5]],
        'Case2': [[0.5217391304347826], [0.0], [0.53125], [0.0], [1.0],
                  [0.0, 0.0], [0.0], [0.5], [0.1], [0.0]],
        'Case3': [[0.475], [0.0], [1.0], [0.4], [0.0], [0.5, 0.0], [1.0],
                  [0.44], [0.15], [0.5]]
    }
    records = get_variants(ad_vcf)
    for c_id, exp in expected.items():
        for i, rec in enumerate(records):
            for j in range(len(rec.alts)):
                result = _get_ad_vaf(rec, c_id, j + 1)
                assert_almost_equals(result, exp[i][j])


def test_fb_vaf_values():
    ''' VAF from AO and RO '''
    expected = {
        'Case1': [[0.4], [0.4666666666666667], [1.0], [0.5], [1.0], [0.0, 0.5],
                  [0.21875], [1.0], [0.0], [0.5]],
        'Case2': [[0.5217391304347826], [0.0], [0.53125], [0.0], [1.0],
                  [0.0, 0.0], [0.0], [0.5], [0.1], [0.0]],
        'Case3': [[0.475], [0.0], [1.0], [0.4], [0.0], [0.5, 0.0], [1.0],
                  [0.44], [0.15], [0.5]]
    }
    records = get_variants(fb_vcf)
    for c_id, exp in expected.items():
        for i, rec in enumerate(records):
            for j in range(len(rec.alts)):
                result = _get_freebayes_vaf(rec, c_id, j + 1)
                assert_almost_equals(result, exp[i][j])


def test_nv_vaf_values():
    ''' VAF from NV and NR '''
    expected = {
        'Case1': [[0.4], [0.4666666666666667], [1.0], [0.5], [1.0], [0.0, 0.5],
                  [0.21875], [1.0], [0.0], [0.5]],
        'Case2': [[0.5217391304347826], [0.0], [0.53125], [0.0], [1.0],
                  [0.0, 0.0], [0.0], [0.5], [0.1], [0.0]],
        'Case3': [[0.475], [0.0], [1.0], [0.4], [0.0], [0.5, 0.0], [1.0],
                  [0.44], [0.15], [0.5]]
    }
    records = get_variants(nv_vcf)
    for c_id, exp in expected.items():
        for i, rec in enumerate(records):
            for j in range(len(rec.alts)):
                result = _get_platypus_vaf(rec, c_id, j + 1)
                assert_almost_equals(result, exp[i][j])


def test_strelka_vaf_values():
    ''' VAF from Strelka somatic VCFs '''
    expected = {
        'Case1': [[0.4], [0.4666666666666667], [1.0], [0.5], [1.0], [0.0, 0.5],
                  [0.21875], [1.0], [0.0], [0.5]],
        'Case2': [[0.5217391304347826], [0.0], [0.53125], [0.0], [1.0],
                  [0.0, 0.0], [0.0], [0.5], [0.1], [0.0]],
        'Case3': [[0.475], [0.0], [1.0], [0.4], [0.0], [0.5, 0.0], [1.0],
                  [0.44], [0.15], [0.5]]
    }
    records = get_variants(strelka_vcf)
    for c_id, exp in expected.items():
        for i, rec in enumerate(records):
            for j in range(len(rec.alts)):
                result = _get_strelka_vaf(rec, c_id, j + 1)
                assert_almost_equals(result, exp[i][j])


def test_svaba_vaf_values():
    ''' VAF from SVABa somatic indels '''
    expected = {
        'Case1': [[0.4], [0.4666666666666667], [1.0], [0.5], [1.0], [0.0],
                  [0.5], [0.21875], [1.0], [0.0], [0.5]],
        'Case2': [[0.5217391304347826], [0.0], [0.53125], [0.0], [1.0],
                  [0.0], [0.0], [0.0], [0.5], [0.1], [0.0]],
        'Case3': [[0.475], [0.0], [1.0], [0.4], [0.0], [0.5], [0.0], [1.0],
                  [0.44], [0.15], [0.5]]
    }
    records = get_variants(svaba_vcf)
    for c_id, exp in expected.items():
        for i, rec in enumerate(records):
            for j in range(len(rec.alts)):
                result = _get_svaba_vaf(rec, c_id, j + 1)
                assert_almost_equals(result, exp[i][j])


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
