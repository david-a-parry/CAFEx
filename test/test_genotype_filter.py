import os
import pysam
from nose.tools import *
from case_control_filter.genotype_filter import FormatFilter

dir_path = os.path.dirname(os.path.realpath(__file__))
ad_vcf = os.path.join(dir_path, 'test_data', 'ad_test.vcf')
nv_vcf = os.path.join(dir_path, 'test_data', 'nv_test.vcf')
fb_vcf = os.path.join(dir_path, 'test_data', 'fb_test.vcf')


def _get_variants(path):
    records = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            records.append(rec)
    return records


def _get_f_filter(expressions, vcf):
    with pysam.VariantFile(vcf) as variants:
        f_filter = FormatFilter(variants, expressions)
    return f_filter


def check_filters(expected, expressions, vcf=ad_vcf):
    records = _get_variants(vcf)
    format_filter = _get_f_filter(expressions, vcf)
    for c_ids, exp in expected.items():
        for i, rec in enumerate(records):
            result = format_filter.filter(rec, c_ids)
            assert_equal(result, exp[i])


def test_gq_single():
    expressions = ['GQ <= 50'.split()]
    expected = {
        ('Case1',): [1, 0, 1, 1, 0, 3, 1, 0, 1, 1],
        ('Case2',): [0, 1, 0, 1, 0, 3, 1, 1, 1, 0],
        ('Case3',): [0, 1, 1, 1, 1, 0, 0, 1, 1, 0]
    }
    check_filters(expected, expressions)


def test_gq_all():
    expressions = ['GQ <= 50 all'.split()]
    expected = {
        ('Case1',): [1, 0, 1, 1, 0, 3, 1, 0, 1, 1],
        ('Case2',): [0, 1, 0, 1, 0, 3, 1, 1, 1, 0],
        ('Case3',): [0, 1, 1, 1, 1, 0, 0, 1, 1, 0],
        ('Case1', 'Case2'): [0, 0, 0, 1, 0, 3, 1, 0, 1, 0],
        ('Case1', 'Case2', 'Case3',): [0, 0, 0, 1, 0, 0, 0, 0, 1, 0]
    }
    check_filters(expected, expressions)


def test_ad_single():
    expressions = ['AD > 3'.split()]
    expected = {
        ('Case1',): [1, 1, 1, 0, 1, 2, 1, 1, 0, 1],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        ('Case3',): [1, 0, 1, 1, 0, 1, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 1, 1, 3, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)


def test_ad_multiple():
    expressions = ['AD > 3 2'.split()]
    expected = {
        ('Case1', 'Case2'): [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        ('Case2', 'Case3'): [1, 0, 1, 0, 0, 0, 0, 1, 0, 0],
        ('Case1', 'Case3'): [1, 0, 1, 0, 0, 0, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 0, 1, 0, 1, 0, 1, 1, 0, 1],
    }
    check_filters(expected, expressions)

# TODO check Number=A (e.g. 'AO' tag from Freebayes)
# TODO check Number=. (e.g. 'NV' tag from Platypus)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
