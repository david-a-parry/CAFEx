import os
import tempfile
from nose.tools import *
from case_control_filter.case_control_filter import main
from .utils import get_variants

dir_path = os.path.dirname(os.path.realpath(__file__))
ad_vcf = os.path.join(dir_path, 'test_data', 'ad_test.vcf')
nv_vcf = os.path.join(dir_path, 'test_data', 'nv_test.vcf')
fb_vcf = os.path.join(dir_path, 'test_data', 'fb_test.vcf')
svaba_vcf = os.path.join(dir_path, 'test_data', 'svaba_test.vcf')
strelka_vcf = os.path.join(dir_path, 'test_data', 'strelka_test.vcf')

ad_records = get_variants(ad_vcf)


def get_tmp_out(suffix='.vcf'):
    f, fname = tempfile.mkstemp(suffix=suffix)
    return fname


def test_missing_samples():
    ''' Raise ValueError on missing samples'''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3', 'Case4'],
                  output=out)
    assert_raises(ValueError, main, ad_vcf, **kwargs)
    kwargs['control'] = ['Control4']
    assert_raises(ValueError, main, ad_vcf, **kwargs)
    if os.path.exists(out):
        os.remove(out)


def test_one_case_genotypes():
    ''' Filter on single case genotypes '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1'], output=out, quiet=True)
    expected_indices = [True, True, False, True, False, True, False, True,
                        False, True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_one_control_genotypes():
    ''' Filter on single control genotypes '''
    raise NotImplementedError()


def test_one_case_control_genotypes():
    ''' Filter on single case/control pair genotypes '''
    raise NotImplementedError()


def test_multiple_case_genotypes():
    ''' Filter on multiple case genotypes '''
    raise NotImplementedError()


def test_multiple_control_genotypes():
    ''' Filter on multiple control genotypes '''
    raise NotImplementedError()


def test_expressions_plus_genotypes():
    ''' Filter cases using both genotypes and expressions '''
    raise NotImplementedError()


def test_expressions_without_genotypes():
    ''' Filter cases on expressions only '''
    raise NotImplementedError()


def test_vaf_filtering():
    ''' Filter on genotypes, VAF and expressions '''
    raise NotImplementedError()


def test_vaf_ratio():
    ''' Filter on case/control VAF ratio '''
    raise NotImplementedError()

