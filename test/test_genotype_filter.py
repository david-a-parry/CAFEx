import os
import pysam
from nose.tools import *
from .utils import get_variants
from cafex.genotype_filter import FormatFilter

dir_path = os.path.dirname(os.path.realpath(__file__))
ad_vcf = os.path.join(dir_path, 'test_data', 'ad_test.vcf')
nv_vcf = os.path.join(dir_path, 'test_data', 'nv_test.vcf')
fb_vcf = os.path.join(dir_path, 'test_data', 'fb_test.vcf')


def _get_f_filter(expressions, vcf):
    with pysam.VariantFile(vcf) as variants:
        f_filter = FormatFilter(variants, expressions)
    return f_filter


def check_filters(expected, expressions, vcf=ad_vcf):
    records = get_variants(vcf)
    format_filter = _get_f_filter(expressions, vcf)
    for c_ids, exp in expected.items():
        for i, rec in enumerate(records):
            result = format_filter.filter(rec, c_ids)
            assert_equal(result, exp[i])


def test_short_expression():
    ''' Raise ValueError if expression too short'''
    assert_raises(ValueError, _get_f_filter, ["AD 2"], ad_vcf)


def test_hanging_expression():
    ''' Raise ValueError if hanging values in expression'''
    assert_raises(ValueError, _get_f_filter, ["AD > 2 DP 2"], ad_vcf)


def test_invalid_min_match():
    ''' Raise ValueError if minimum matching parameter is not a number > 0'''
    assert_raises(ValueError, _get_f_filter, ["AD > 2 foo"], ad_vcf)
    assert_raises(ValueError, _get_f_filter, ["AD > 2 0"], ad_vcf)
    assert_raises(ValueError, _get_f_filter, ["AD > 2 -1"], ad_vcf)


def test_invalid_operator():
    ''' Raise ValueError if invalid operator used'''
    assert_raises(ValueError, _get_f_filter, ["AD ^ 2"], ad_vcf)


def test_invalid_logical_operator():
    ''' Raise ValueError if invalid logical operator used'''
    assert_raises(ValueError, _get_f_filter, ["AD > 2 xor DP < 20"], ad_vcf)


def test_invalid_format_field():
    ''' Raise ValueError if invalid FORMAT field used'''
    assert_raises(ValueError, _get_f_filter, ["XX > 2"], ad_vcf)


def test_invalid_value():
    ''' Raise ValueError if invalid value type used'''
    assert_raises(ValueError, _get_f_filter, ["AD != foo"], ad_vcf)


def test_invalid_sum():
    ''' Raise ValueError if attempt to sum non-subscriptable field'''
    assert_raises(ValueError, _get_f_filter, ["sum(GQ) > 10"], ad_vcf)


def test_invalid_subscript():
    ''' Raise ValueError if attempt to subscript non-subscriptable field'''
    assert_raises(ValueError, _get_f_filter, ["GQ[1] > 10"], ad_vcf)
    assert_raises(ValueError, _get_f_filter, ["GQ[0] > 10"], ad_vcf)


def test_gq_single():
    ''' Number=1 '''
    expressions = ['GQ <= 50']
    expected = {
        ('Case1',): [1, 0, 1, 1, 0, 3, 1, 0, 1, 1],
        ('Case2',): [0, 1, 0, 1, 0, 3, 1, 1, 1, 0],
        ('Case3',): [0, 1, 1, 1, 1, 0, 0, 1, 1, 0]
    }
    check_filters(expected, expressions)


def test_ao_single():
    ''' Number=A '''
    expressions = ['AO > 3']
    expected = {
        ('Case1',): [1, 1, 1, 0, 1, 2, 1, 1, 0, 1],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        ('Case3',): [1, 0, 1, 1, 0, 1, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 1, 1, 3, 1, 1, 0, 1],
    }
    check_filters(expected, expressions, vcf=fb_vcf)


def test_nv_single():
    ''' Number=. '''
    expressions = ['NV > 3']
    expected = {
        ('Case1',): [1, 1, 1, 0, 1, 3, 1, 1, 0, 1],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        ('Case3',): [1, 0, 1, 1, 0, 3, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 1, 1, 3, 1, 1, 0, 1],
    }
    check_filters(expected, expressions, vcf=nv_vcf)


def test_ad_single():
    ''' Number=R '''
    expressions = ['AD > 3']
    expected = {
        ('Case1',): [1, 1, 1, 0, 1, 2, 1, 1, 0, 1],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        ('Case3',): [1, 0, 1, 1, 0, 1, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 1, 1, 3, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)


def test_gq_all():
    ''' Number=1 all samples matching '''
    expressions = ['GQ <= 50 all']
    expected = {
        ('Case1',): [1, 0, 1, 1, 0, 3, 1, 0, 1, 1],
        ('Case2',): [0, 1, 0, 1, 0, 3, 1, 1, 1, 0],
        ('Case3',): [0, 1, 1, 1, 1, 0, 0, 1, 1, 0],
        ('Case1', 'Case2'): [0, 0, 0, 1, 0, 3, 1, 0, 1, 0],
        ('Case1', 'Case2', 'Case3',): [0, 0, 0, 1, 0, 0, 0, 0, 1, 0]
    }
    check_filters(expected, expressions)


def test_ad_multiple():
    ''' Number=R with multiple samples matching '''
    expressions = ['AD > 3 2']
    expected = {
        ('Case1', 'Case2'): [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        ('Case2', 'Case3'): [1, 0, 1, 0, 0, 0, 0, 1, 0, 0],
        ('Case1', 'Case3'): [1, 0, 1, 0, 0, 0, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 0, 1, 0, 1, 0, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)


def test_two_expressions():
    ''' Multiple expressions '''
    expressions = ['AD > 3', 'GQ > 50']
    expected = {
        ('Case1',): [0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        ('Case3',): [1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 0, 1, 3, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)


def test_gq_or_ad():
    ''' Logical or '''
    expressions = ['AD > 3 or GQ > 50']
    expected = {
        ('Case1',): [1, 1, 1, 0, 1, 2, 1, 1, 0, 1],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 1, 0, 1],
        ('Case3',): [1, 0, 1, 1, 0, 3, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 1, 1, 3, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)


def test_gq_and_ad():
    ''' Logical and '''
    expressions = ['AD > 3 and GQ > 50']
    expected = {
        ('Case1',): [0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        ('Case3',): [1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 0, 1, 1, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)


def test_gq_and_ad_or_dp():
    ''' Logical and plus or '''
    expressions = ['AD > 3 and GQ > 50 or DP > 20']
    expected = {
        ('Case1',): [0, 1, 0, 0, 1, 3, 1, 1, 0, 0],
        ('Case2',): [1, 0, 1, 1, 1, 3, 0, 0, 1, 1],
        ('Case3',): [1, 0, 1, 0, 0, 3, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 1, 1, 3, 1, 1, 1, 1],
    }
    check_filters(expected, expressions)


def test_gq_and_ad_and_dp():
    ''' Logical and times 2 '''
    expressions = ['AD > 3 and GQ > 50 and DP > 20']
    expected = {
        ('Case1',): [0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        ('Case3',): [1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 0, 1, 0, 1, 1, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)
    expressions = ['AD > 3 and GQ > 50 and DP > 20 2']
    expected = {
        ('Case1', 'Case2', 'Case3'): [1, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    }
    check_filters(expected, expressions)


def test_sum_ad():
    ''' Sum value comparison '''
    expressions = ['sum(AD) == 20']
    expected = {
        ('Case1',): [1, 0, 1, 0, 0, 0, 0, 0, 1, 1],
        ('Case2',): [0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
        ('Case3',): [0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
        ('Case1', 'Case2', 'Case3'): [1, 0, 1, 0, 1, 0, 1, 1, 1, 1],
    }
    check_filters(expected, expressions)
    expressions = ['sum(AD) > 20 2']
    expected = {
        ('Case1', 'Case2', 'Case3'): [1, 0, 1, 0, 1, 3, 1, 1, 0, 1]
    }
    check_filters(expected, expressions)


def test_subscript_ad():
    ''' Subscripted value comparison '''
    expressions = ['AD[1] == 30']
    expected = {
        ('Case1',): [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        ('Case2',): [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        ('Case3',): [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    }
    check_filters(expected, expressions)
    expressions = ['AD[1] >= 20 2']
    expected = {
        ('Case1', 'Case2', 'Case3'): [0, 0, 1, 0, 1, 0, 0, 0, 0, 0]
    }
    check_filters(expected, expressions)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
