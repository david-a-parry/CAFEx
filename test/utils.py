import pysam


def get_variants(path):
    records = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            records.append(rec)
    return records
