#!/usr/bin/env python3
import pysam
import os
from random import randint

dir_path = os.path.dirname(os.path.realpath(__file__))
ad_vcf = os.path.join(dir_path, 'test_data', 'ad_test.vcf')


def ad_to_nv(record):
    new_record = record.copy()
    for sample in record.samples:
        ads = record.samples[sample]['AD']
        nv = ads[1:]
        nr = [ads[0] + sum(nv)] * (len(ads) - 1)
        new_record.samples[sample]['NR'] = nr
        new_record.samples[sample]['NV'] = nv
    del new_record.format['AD']
    del new_record.format['DP']
    return new_record


def ad_to_fb(record):
    new_record = record.copy()
    for sample in record.samples:
        ads = record.samples[sample]['AD']
        ro = ads[0]
        ao = ads[1:]
        new_record.samples[sample]['RO'] = ro
        new_record.samples[sample]['AO'] = ao
    del new_record.format['AD']
    return new_record


def ad_to_strelka(record):
    new_record = record.copy()
    for sample in record.samples:
        ads = record.samples[sample]['AD']
        if len(record.ref) == len(record.alts[0]) == 1:
            for nt in 'ACGT':
                if nt == record.ref:
                    n = ads[0]
                elif nt in record.alts:
                    i = record.alleles.index(nt)
                    n = ads[i]
                else:
                    n = 0
                new_record.samples[sample][nt + 'U'] = (n, n + randint(0, 11))
        else:
            new_record.samples[sample]['TAR'] = (ads[0], ads[0] + randint(0, 9))
            new_record.samples[sample]['TIR'] = (ads[1], ads[1] + randint(0, 9))
    del new_record.format['AD']
    return new_record


def ad_to_svaba(record):
    new_records = []
    for i in range(1, len(record.alleles)):
        new_rec = record.copy()
        new_rec.alleles = (record.alleles[0], record.alleles[i])
        for sample in record.samples:
            ads = record.samples[sample]['AD']
            # hack to get around pysam complaining that it expects two values
            new_rec.samples[sample]['AD'] = (ads[i], None)
            new_rec.samples[sample]['GT'] = tuple(1 if x == i else 0 for x in
                                                  record.samples[sample]['GT'])
        new_rec = str(new_rec).replace(',.:', ':')
        new_records.append(new_rec)
    return new_records


conversions = {
    'nv': {'vcf': os.path.join(dir_path, 'test_data', 'nv_test.vcf'),
           'method': ad_to_nv},
    'fb': {'vcf': os.path.join(dir_path, 'test_data', 'fb_test.vcf'),
           'method': ad_to_fb},
    'strelka': {'vcf': os.path.join(dir_path, 'test_data', 'strelka_test.vcf'),
                'method': ad_to_strelka},
    'svaba': {'vcf': os.path.join(dir_path, 'test_data', 'svaba_test.vcf'),
              'method': ad_to_svaba}
}


def get_outputs(input_vcf):
    output_vcfs = dict()
    for k, d in conversions.items():
        vcf = d['vcf']
        header = vcf.replace('.vcf', '.header.vcf')
        with pysam.VariantFile(header) as hd:
            output_vcfs[k] = open(vcf, mode='wt')
            output_vcfs[k].write(str(hd.header))
            if k == 'nv':
                add = ('NR', 'NV')
            elif k == 'fb':
                add = ('RO', 'AO')
            elif k == 'strelka':
                add = ('AU', 'CU', 'GU', 'TU', 'TAR', 'TIR')
            elif k == 'svaba':
                add = ()
            else:
                raise ValueError("Don't recognise key '{}'".format(k))
            for f in add:
                fmt = hd.header.formats[f]
                input_vcf.header.formats.add(f,
                                             fmt.number,
                                             fmt.type,
                                             fmt.description)
    return output_vcfs


def main():
    with pysam.VariantFile(ad_vcf) as ad_vars:
        output_vcfs = get_outputs(ad_vars)
        for record in ad_vars:
            for k, out in output_vcfs.items():
                converted = conversions[k]['method'](record)
                if isinstance(converted, list):
                    for conv in converted:
                        out.write(str(conv))
                else:
                    out.write(str(converted))
    for out in output_vcfs.values():
        out.close()


if __name__ == '__main__':
    main()
