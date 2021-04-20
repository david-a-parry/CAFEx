import logging
import pysam
import sys
import time
from case_control_filter.bit_utils import set_first_bits
from case_control_filter.vaf import get_vaf_method
from case_control_filter.genotype_filter import FormatFilter

PROG_NAME = "case_control_filter"
logger = logging.getLogger(PROG_NAME)
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)


def check_samples(vcf, samples):
    missing = [x for x in samples if x not in vcf.header.samples]
    if missing:
        raise ValueError("The following specified samples were not found in " +
                         "the VCF: " + ",".join(missing))


def get_format_filter(vcf, expressions):
    if not expressions:
        return None
    return FormatFilter(vcf=vcf, expressions=expressions)


def get_max_vafs(record, samples, vaf_func):
    vafs = []
    for i in range(1, len(record.alleles)):
        vafs.append(max(vaf_func(record, x, i) for x in samples))
    return vafs


def add_info_tag(vcf, tag):
    if tag in vcf.header.info:
        if vcf.header.info[tag].number != 1:
            raise ValueError("INFO tag '{}' already exists ".format(tag) +
                             "but has Number '{}' - will not overwite."
                             .format(vcf.header.info[tag].number))
        logger.warn("Overwriting pre-existing '{}' INFO field".format(tag))
    vcf.header.info.add(tag, '1', 'Integer',
                        'Bitwise flag indicating which values passed filter ' +
                        'expressions from ' + PROG_NAME)


def main(vcf, case=[], control=[], output=None, ignore_genotypes=False,
         case_expressions=[], control_expressions=[], min_case_vaf=None,
         max_control_vaf=None, vaf_ratio=None, info_tag=None,
         progress_interval=100_000, quiet=False, debug=False):
    if quiet:
        logger.setLevel(logging.WARN)
    elif debug:
        logger.setLevel(logging.DEBUG)
    output = '-' if output is None else output
    with pysam.VariantFile(vcf) as variants:
        check_samples(variants, case + control)
        vaf_calculation = None
        if min_case_vaf or max_control_vaf or vaf_ratio:
            vaf_calculation = get_vaf_method(variants)
        case_filter = get_format_filter(variants, case_expressions)
        control_filter = get_format_filter(variants, control_expressions)
        if info_tag:
            add_info_tag(variants, info_tag)
        out = pysam.VariantFile(output, 'w', header=variants.header)
        out.header.add_meta(key=PROG_NAME,
                            value=str.join(" ", sys.argv) + "; Date=" +
                            time.strftime("%Y-%m-%d %H:%M"))
        read, written = 0, 0
        for record in variants:
            if progress_interval and read:
                if read % progress_interval == 0:
                    logger.info(
                        "{:,} variants processed, {:,} written. At {}:{}"
                        .format(read, written, record.chrom, record.pos))
            read += 1
            n_alts = len(record.alts)
            filter_flag = set_first_bits(n_alts)  # bitwise flag per ALT
            for alt in range(n_alts):
                allele = alt + 1
                if not ignore_genotypes:
                    if all(allele not in record.samples[x]['GT'] for x in
                           case):
                        filter_flag &= ~(1 << alt)  # unset bit for ALT allele
                    if any(allele in record.samples[x]['GT'] for x in control):
                        filter_flag &= ~(1 << alt)
            # are we done already?
            if not filter_flag:
                continue
            if control_filter is not None:
                filter_flag &= control_filter.filter(record, control)
            if not filter_flag:
                continue
            if case_filter is not None:
                filter_flag &= case_filter.filter(record, case)
            if not filter_flag:
                continue
            ca_vafs = None
            co_vafs = None
            if max_control_vaf or vaf_ratio:
                co_vafs = get_max_vafs(record, control, vaf_calculation)
                if max_control_vaf:
                    for i in range(n_alts):
                        filter_flag &= ~((co_vafs[i] > max_control_vaf) << i)
                    if not filter_flag:
                        continue
            if min_case_vaf or vaf_ratio:
                ca_vafs = get_max_vafs(record, case, vaf_calculation)
                if min_case_vaf:
                    for i in range(n_alts):
                        filter_flag &= ~((ca_vafs[i] < min_case_vaf) << i)
                    if not filter_flag:
                        continue
            if vaf_ratio:
                for i in range(n_alts):
                    if co_vafs[i] > 0.0:
                        filter_flag &= ~((ca_vafs[i]/co_vafs[i] < vaf_ratio)
                                         << i)
                    elif ca_vafs[i] == 0.0:
                        filter_flag &= ~(1 << i)  # unset if case VAF is zero
            if filter_flag:
                if info_tag:
                    record.info[info_tag] = filter_flag
                out.write(record)
                written += 1
    logger.info("Finished processing {:,} variants. ".format(read) +
                "{:,} written, {:,} filtered.".format(written, read - written))
    out.close()
