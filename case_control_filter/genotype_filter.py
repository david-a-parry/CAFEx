import operator
from .bit_utils import first_n_bits_set, set_first_bits, flag_consensus

ops = {
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
    "=": operator.eq,
    "!": operator.ne,
}


class FormatFilter(object):
    '''
        A class for filtering on given FORMAT fields in a VCF
    '''

    def __init__(self, vcf, filters):
        '''
            Args:
                vcf:    VariantFile object from pysam.

                filters:
                        iterable of tuples of field names, operands and values
                        for filtering. Optional fourth parameter indicates how
                        many samples parameter has to match (either an integer
                        or 'all' are acceptable). By default only one sample
                        passed to the 'filter' method needs match an
                        expression.

        '''
        self.vcf = vcf
        self.metadata = self.vcf.header.formats
        self.fields = set()
        self.filters = self._parse_filters(filters)

    def _parse_filters(self, filters):
        '''
            Args:
                filters:
                    iterable of tuples of field names, operands and values for
                    filtering. Optional fourth parameter indicates how many
                    samples parameter has to match (either an integer or 'all'
                    are acceptable). By default only one sample passed to the
                    'filter' method needs match an expression.
        '''
        f = []
        for field, operand, value, *extra in filters:
            exp = " ".join([field, operand, value] + extra)
            min_matching = 1
            if extra:
                if extra[0].lower() == 'all':
                    min_matching = -1
                else:
                    try:
                        min_matching = int(extra[0])
                    except ValueError:
                        raise ValueError("Invalid minimum matching parameter" +
                                         "'{}'".format(extra[0]) + "passed " +
                                         "to expression '{}'. ".format(exp) +
                                         "Must either be numeric or 'all'.")
                    if min_matching < 1:
                        raise ValueError("Invalid minimum matching parameter" +
                                         " '{}'".format(extra[0]) + "passed " +
                                         "to expression '{}'. ".format(exp) +
                                         "Must be > 0.")
            try:
                op = ops[operand]
            except KeyError:
                raise ValueError("Unrecognised operand '{}' ".format(operand) +
                                "in filter expression '{} {} {}'".format(
                                field, operand, value))
            if field not in self.metadata:
                raise ValueError("FORMAT field '{}' not in VCF ".format(field) +
                                 "header - can not be used for FORMAT field " +
                                 "filtering.")
            ftype = self.metadata[field].type
            coerc = None
            if ftype == 'Integer':
                coerc = int
            elif ftype == 'Float':
                coerc = float
            if coerc is not None:
                try:
                    value = coerc(value)
                except ValueError:
                    raise ValueError("Filter value for FORMAT field " +
                                     "'{}' could not be ".format(field) +
                                     "converted to {}, but".format(ftype) +
                                     "but field Type is {}".format(ftype) +
                                     "in VCF header.")
            num = self.metadata[field].number
            f.append((field, op, value, num, min_matching))
            self.fields.add(field)
        return f

    def _check_sample(self, record, smpl, field, op, value, number, n_alts):
        '''
            Return bitwise flag indicating for each alt allele whether sample
            call meets expression criteria.
        '''
        if field not in record.samples[smpl]:
            return 0
        flag = 0  # set nth bit to 1 if nth ALT matches expression
        annot = record.samples[smpl][field]
        if number == 'A' or number == 'R':
            if number == 'R':
                alt_vals = annot[1:]
            else:
                alt_vals = annot
            for i in range(n_alts):
                if alt_vals[i] is not None and op(alt_vals[i], value):
                    flag |= 1 << i  # set bit for ALT allele
        else:
            if number == 1:
                if annot is not None and op(annot, value):
                    return set_first_bits(n_alts)
            else:
                for x in annot:
                    if x is not None and op(x, value):
                        # set all bits if ANY value matches
                        return set_first_bits(n_alts)
        return flag

    def filter(self, record, samples):
        '''
            Read VcfRecord and returns a bitwsise flag indicating
            whether each ALT allele meets parameters from self.filters.

            record:    VcfRecord to assess according to self.filters.

            samples:   Iterable of samples to assess.

        '''
        n_alts = len(record.alts)
        # array of bitwise flags, one flag per filter, bits set if allele meets
        # filter criteria
        flags = []
        for field, op, value, number, min_smpls in self.filters:
            alt_f = 0  # per-alt flag for this filter expression
            alt_counts = None
            if min_smpls < 1:
                min_smpls = len(samples)
            if min_smpls > 1:
                alt_counts = [0] * n_alts
            for smp in samples:
                flt = self._check_sample(record, smp, field, op, value, number,
                                         n_alts)
                if min_smpls == 1:
                    if first_n_bits_set(flt, n_alts):
                        alt_f = flt  # bail out early if all alleles pass
                        break
                    alt_f |= flt
                else:
                    for i in range(n_alts):
                        alt_counts[i] += flt >> i & 1
            if alt_counts is not None:
                for i in range(n_alts):
                    alt_f |= (alt_counts[i] >= min_smpls) << i
            flags.append(alt_f)
        # report whether ALL filter expressions were matched
        return flag_consensus(flags)
