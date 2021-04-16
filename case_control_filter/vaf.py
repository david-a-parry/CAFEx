def _vaf_ad_dp(ad, dp):
    if dp > 0.0:
        return ad/dp
    return 0.0


def _get_ad_vaf(record, sample, allele):
    ad = record.samples[sample]['AD']
    dp = sum(ad)
    return _vaf_ad_dp(ad[allele], dp)


def _get_svaba_vaf(record, sample, allele):
    # only ever 1 value for AD in SvABA
    ad = record.samples[sample]['AD']
    dp = record.samples[sample]['DP']
    return _vaf_ad_dp(ad, dp)


def _get_strelka_snv_vaf(record, sample, allele):
    ref_k = record.ref + 'U'
    alt_k = record.alleles[allele] + 'U'
    ad = record.samples[sample][alt_k][0]
    dp = record.samples[sample][ref_k][0] + ad
    return _vaf_ad_dp(ad, dp)


def _get_strelka_indel_vaf(record, sample, allele):
    ad = record.samples[sample]['TIR'][0]
    dp = record.samples[sample]['TAR'][0] + ad
    return _vaf_ad_dp(ad, dp)


def _get_strelka_vaf(record, sample, allele):
    if 'AU' in record.format:
        return _get_strelka_snv_vaf(record, sample, allele)
    return _get_strelka_indel_vaf(record, sample, allele)


def _get_platypus_vaf(record, sample, allele):
    ad = record.samples[sample]['NV'][allele - 1]
    dp = record.samples[sample]['NR'][allele - 1]
    return _vaf_ad_dp(ad, dp)


def _get_freebayes_vaf(record, sample, allele):
    ad = record.samples[sample]['AO'][allele - 1]
    dp = record.samples[sample]['RO'][0] + ad
    return _vaf_ad_dp(ad, dp)


def get_vaf_method(vcf):
    '''
    Scan VCF header to determine which method to use to calculate VAF. Returns
    a function to calculate VAF for given VcfRecord, sample and allele index.

    Will use AD field if found but non-standard fields from Strelka, Platypus
    and Freebayes are also supported.

    Example:
            vcf = pysam.VariantFile('input.bcf')
            vaf_calc = get_vaf_method(vcf)
            for record in vcf:
                sample_vaf = vaf_calc(record, 'Sample1', 1)

    '''
    if 'AD' in vcf.header.formats:
        if vcf.header.formats['AD'].number == 1:  # SvABA
            return _get_svaba_vaf
        return _get_ad_vaf
    if 'AU' in vcf.header.formats or 'TAR' in vcf.header.formats:  # Strelka
        # usually Strelka SNVs and Indels will be in separate VCFs
        if 'AU' in vcf.header.formats and 'TAR' in vcf.header.formats:
            # presumably combined VCF
            return _get_strelka_vaf
        elif 'AU' in vcf.header.formats:
            return _get_strelka_snv_vaf
        return _get_strelka_indel_vaf
    if 'NV' in vcf.header.formats and 'NR' in vcf.header.formats:
        return _get_platypus_vaf
    if 'AO' in vcf.header.formats and 'RO' in vcf.header.formats:
        return _get_freebayes_vaf
    raise ValueError("Could not identify any supported allele depth fields " +
                     "in input VCF header.")
