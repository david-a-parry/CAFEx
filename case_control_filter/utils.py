def all_bits_set(x, nbits):
    ''' Return True if first nbits bits are set. '''
    return all(x & 1 << y for y in range(nbits))


def set_bits_in_range(nbits):
    ''' Set first n bits to 1 '''
    x = 0
    for i in range(nbits):
        x |= (1 << i)
    return x


def flip_bits(x, nbits):
    '''Flip the first nbits bits of x'''
    return x ^ set_bits_in_range(nbits)


def highest_set_bit(x):
    ''' Return index of highest set bit in x'''
    if x == 0:
        return 0
    hsb = 0
    x = int(x / 2)
    while x > 0:
        x = int(x / 2)
        hsb += 1
    return hsb


def flag_consensus(flags):
    ''' For list of flags return bitwise & for all. '''
    f = flags[0]
    for x in flags[1:]:
        f &= x
    return f


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
    dp = record.samples[sample]['NR'][0] + ad
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

