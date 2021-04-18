import os
import pysam
from nose.tools import *
from case_control_filter.genotype_filter import FormatFilter

dir_path = os.path.dirname(os.path.realpath(__file__))
input_vcf = os.path.join(dir_path, 'test_data', 'ad_test.vcf')


def _get_variants(path):
    records = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            records.append(rec)
    return records


def _get_f_filter(expressions):
    with pysam.VariantFile(input_vcf) as vcf:
        f_filter = FormatFilter(vcf, expressions)
    return f_filter


def test_ad_single():
    expressions = ['AD > 3'.split()]
    format_filter = _get_f_filter(expressions)
    records = _get_variants(input_vcf)
    expected = {
        ('Case1',): [1, 1, 1, 0, 1, 2, 1, 1, 0, 1],
        ('Case2',): [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        ('Case3',): [1, 0, 1, 1, 0, 1, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 1, 1, 1, 1, 3, 1, 1, 0, 1],
    }
    for c_ids, exp in expected.items():
        for i, rec in enumerate(records):
            result = format_filter.filter(rec, c_ids)
            assert_equal(result, exp[i])


def test_ad_multiple():
    expressions = ['AD > 3 2'.split()]
    format_filter = _get_f_filter(expressions)
    records = _get_variants(input_vcf)
    expected = {
        ('Case1', 'Case3'): [1, 0, 1, 0, 0, 0, 1, 1, 0, 1],
        ('Case1', 'Case2', 'Case3'): [1, 0, 1, 0, 1, 0, 1, 1, 0, 1],
    }
    for c_ids, exp in expected.items():
        for i, rec in enumerate(records):
            result = format_filter.filter(rec, c_ids)
            assert_equal(result, exp[i])


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
