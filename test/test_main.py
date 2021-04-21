import os
import tempfile
from nose.tools import *
from cafex.case_control_filter import main
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


def test_no_samples():
    ''' Raise ValueError if no samples provided '''
    out = get_tmp_out()
    kwargs = dict(case=[], control=[], output=out)
    assert_raises(ValueError, main, ad_vcf, **kwargs)
    if os.path.exists(out):
        os.remove(out)


def test_missing_samples():
    ''' Raise ValueError on missing samples'''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3', 'Case4'],
                  output=out)
    assert_raises(ValueError, main, ad_vcf, **kwargs)
    kwargs = dict(control=['Control4'],
                  output=out)
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
    out = get_tmp_out()
    kwargs = dict(control=['Control1'], output=out, quiet=True)
    expected_indices = [False, True, True, True, True, True, True, True, True,
                        True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_one_case_control_genotypes():
    ''' Filter on single case/control pair genotypes '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1'],
                  control=['Control1'],
                  output=out,
                  quiet=True)
    expected_indices = [False, True, False, True, False, True, False, True,
                        False, True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_multiple_case_genotypes():
    ''' Filter on multiple case genotypes '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  output=out,
                  quiet=True)
    expected_indices = [True, True, True, True, True, True, False, True, True,
                        True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_multiple_control_genotypes():
    ''' Filter on multiple control genotypes '''
    out = get_tmp_out()
    kwargs = dict(control=['Control1', 'Control2', 'Control3'],
                  output=out,
                  quiet=True)
    expected_indices = [False, True, True, False, True, True, True, True, True,
                        True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_multiple_case_control_genotypes():
    ''' Filter on multiple case/control genotypes '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  control=['Control1', 'Control2', 'Control3'],
                  output=out,
                  quiet=True)
    expected_indices = [False, True, True, False, True, True, False, True,
                        True, True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_expression_plus_genotypes():
    ''' Filter cases using both genotypes and an expression '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  case_expressions=["AD > 5"],
                  output=out,
                  quiet=True)
    expected_indices = [True, True, True, False, True, True, False, True,
                        False, True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_multiexpressions_plus_genotypes():
    ''' Filter cases using both genotypes and multiple expressions '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', ],
                  case_expressions=["DP > 25", "GQ > 50"],
                  output=out,
                  quiet=True)
    expected_indices = [False, False, True, False, True, False, False, True,
                        False, True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_expressions_without_genotypes():
    ''' Filter cases on expressions only '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  case_expressions=["AD > 5"],
                  ignore_genotypes=True,
                  output=out,
                  quiet=True)
    expected_indices = [True, True, True, False, True, True, True, True, False,
                        True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_expressions_with_subscripted_field():
    ''' Filter cases using a subscripted field '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  case_expressions=["AD > 5"],
                  ignore_genotypes=True,
                  output=out,
                  quiet=True)
    expected_indices = [True, True, True, False, True, True, True, True, False,
                        True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_expressions_with_summed_values():
    ''' Filter cases using a summed values '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  case_expressions=["sum(AD) >= 30"],
                  ignore_genotypes=True,
                  output=out,
                  quiet=True)
    expected_indices = [True, False, True, False, True, True, True, True, True,
                        True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_case_control_expression_plus_genotypes():
    '''
    Filter with cases and controls using both genotypes and an expression
    '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  case_expressions=["AD > 5"],
                  control=['Control1', 'Control2', 'Control3'],
                  control_expressions=["AD < 3 all", "DP >= 30 1"],
                  output=out,
                  quiet=True)
    expected_indices = [False, False, False, False, True, True, False, True,
                        False, False]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_min_case_vaf_filtering():
    ''' Filter on genotypes, expressions and minimum case VAF '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  case_expressions=["DP > 20"],
                  control=['Control1', 'Control2', 'Control3'],
                  min_case_vaf=0.2,
                  output=out,
                  quiet=True)
    expected_indices = [False, False, False, False, True, True, False, True,
                        False, True]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)


def test_max_control_vaf_filtering():
    ''' Filter on genotypes, expressions and maximum control VAF '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  case_expressions=["DP > 20"],
                  control=['Control1', 'Control2', 'Control3'],
                  max_control_vaf=0.05,
                  output=out,
                  quiet=True)
    expected_indices = [False, False, False, False, True, True, False, True,
                        True, False]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)


def test_case_control_vaf():
    ''' Filter on minimum case and maximum control VAF '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2'],
                  control=['Control1', 'Control2'],
                  min_case_vaf=0.2,
                  max_control_vaf=0.05,
                  ignore_genotypes=True,
                  output=out,
                  quiet=True)
    expected_indices = [False, False, False, True, True, False, False, True,
                        False, False]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)


def test_vaf_ratio():
    ''' Filter on case/control VAF ratio '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1'],
                  control=['Control1'],
                  vaf_ratio=10.0,
                  ignore_genotypes=True,
                  output=out,
                  quiet=True)
    expected_indices = [False, False, False, True, True, True, True, True,
                        False, False]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_vaf_ratio_multi():
    ''' Filter on multiple case/control VAF ratio '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2'],
                  control=['Control1', 'Control2'],
                  vaf_ratio=10.0,
                  ignore_genotypes=True,
                  output=out,
                  quiet=True)
    expected_indices = [False, False, False, True, True, True, False, True,
                        True, False]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    assert_equal(results, expected_records)
    if os.path.exists(out):
        os.remove(out)


def test_add_info_tag():
    ''' Test filtering and adding bitwise filter tag '''
    out = get_tmp_out()
    kwargs = dict(case=['Case1', 'Case2', 'Case3'],
                  control=['Control1', 'Control2', 'Control3'],
                  output=out,
                  info_tag="TEST_TAG",
                  quiet=True)
    expected_indices = [False, True, True, False, True, True, False, True,
                        True, True]
    expected_flags = [1, 1, 1, 3, 1, 1, 1]
    expected_records = [x for x, y in zip(ad_records, expected_indices) if y]
    main(ad_vcf, **kwargs)
    results = get_variants(out)
    for res, exp, flag in zip(results, expected_records, expected_flags):
        for attr in ('pos', 'ref', 'alts'):
            assert_equal(getattr(res, attr), getattr(exp, attr))
        assert_equal(res.info['TEST_TAG'], flag)
    if os.path.exists(out):
        os.remove(out)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)
