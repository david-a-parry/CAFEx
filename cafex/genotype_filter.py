import operator
import re
from collections import namedtuple
from .bit_utils import first_n_bits_set, set_first_bits, flag_consensus

Expression = namedtuple("Expression",
                        "field operator value number sum subscript")

_subscript_re = re.compile(r'(\w+)\[(\d+)\]')
_sum_re = re.compile(r'sum\((\w+)\)')

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

logical = {  # we use bitwise operators because we will be comparing flags
    "and": operator.iand,
    "&&": operator.iand,
    "or": operator.ior,
    "||": operator.ior,
}


class FilterExpression(object):
    ''' A class for holding expression logic for testing genotype values'''

    def __init__(self, expression, vcf):
        '''
        Args:
            expression:
                Expression string that must consist of FORMAT field names,
                operators and values for filtering. These triplets can be
                optionally joined by 'and' or 'or' logical parameters if more
                than one evaluation is to be performed. Optional final
                parameter indicates how many samples must match the expression
                (either an integer or 'all' are acceptable, defaults to 1).
        '''
        self.expressions = []
        self.logical_ops = []
        self.min_samples = 1
        self.metadata = vcf.header.formats
        self._parse_expressions(expression)

    def check_sample(self, record, smpl, n_alts):
        '''
            Return bitwise flag indicating for each alt allele whether sample
            call meets expression criteria.
        '''
        flags = []
        for exp in self.expressions:
            if exp.field not in record.samples[smpl]:
                flags.append(0)
                continue
            flg = 0  # set nth bit to 1 if nth ALT matches expression
            annot = record.samples[smpl][exp.field]
            val = None
            if exp.number == 1:
                val = annot
            elif exp.sum:
                val = sum(filter(None, annot))
            elif exp.subscript is not None:
                val = annot[exp.subscript]
            if val is not None and exp.operator(val, exp.value):
                flg = set_first_bits(n_alts)
            elif exp.number == 1:  # number is one and value is None
                pass
            elif exp.number == 'A' or exp.number == 'R':
                if exp.number == 'R':
                    alt_vals = annot[1:]
                else:
                    alt_vals = annot
                for i in range(n_alts):
                    if alt_vals[i] is not None and exp.operator(alt_vals[i],
                                                                exp.value):
                        flg |= 1 << i  # set bit for ALT allele
            else:
                for x in annot:
                    if x is not None and exp.operator(x, exp.value):
                        # set all bits if ANY value matches
                        flg = set_first_bits(n_alts)
                        break
            flags.append(flg)
        if len(flags) == 1:
            return flags[0]
        fcmp = self.logical_ops[0](flags[0], flags[1])
        for i in range(1, len(self.logical_ops)):
            fcmp = self.logical_ops[i](fcmp, flags[i + 1])
        return fcmp

    def _parse_expressions(self, expression):
        split = expression.split()
        if len(split) < 3:
            raise ValueError("At least three whitespace-separated values are" +
                             " required for expressions.")
        i = 0
        while i < len(split) - 2:
            field, oprtr, value = split[i:i+3]
            try:
                op = ops[oprtr]
            except KeyError:
                raise ValueError("Unrecognised operator '{}' ".format(oprtr) +
                                 "in filter expression '{}'".format(expression)
                                 )
            do_sum = False
            subscript = None
            sub_match = _subscript_re.match(field)
            sum_match = _sum_re.match(field)
            if sub_match is not None:
                field = sub_match.group(1)
                subscript = int(sub_match.group(2))
            elif sum_match is not None:
                field = sum_match.group(1)
                do_sum = True
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
            if num == 1:
                if do_sum:
                    raise ValueError("Cannot use 'sum' in expression '" +
                                     expression + "'. Only expect one value" +
                                     "(Number=1) for field '{}'.".format(field)
                                     )
                if subscript is not None:
                    raise ValueError("Cannot use 'subscript' in expression '" +
                                     expression + "'. Only expect one value" +
                                     "(Number=1) for field '{}'.".format(field)
                                     )
                elif num not in ('A', 'R', '.') and subscript:
                    if subscript > num - 1:
                        raise ValueError("Subscript '[{}]'".format(subscript) +
                                         " out of range for field '{}' in "
                                         .format(field) + " expression '" +
                                         expression + "'. Expect {} values."
                                         .format(num))
            self.expressions.append(Expression(field, op, value, num, do_sum,
                                               subscript))
            if len(split) >= i + 7:  # logical operator + another expression
                logic = split[i+3]
                if logic not in logical:
                    raise ValueError("Could not parse expression '{}'. "
                                     .format(expression) + "Invalid argument" +
                                     "'{}' at position {}. ".format(split[i+3],
                                                                    i + 4) +
                                     "Expected logical operator (and/or).")
                self.logical_ops.append(logical[logic])
                i += 4
            elif i + 3 <= len(split) <= i + 4:  # end of expr or min match
                i += 3
            else:
                raise ValueError("Hanging values at end of expression '{}'."
                                 .format(expression))
        if i < len(split):
            if i != len(split) - 1:
                raise ValueError("Hanging values at end of expression '{}'."
                                 .format(expression))
            if split[i].lower() == 'all':
                min_matching = -1
            else:
                try:
                    min_matching = int(split[i])
                except ValueError:
                    raise ValueError("Invalid minimum matching parameter" +
                                     "'{}'".format(split[i]) + "passed " +
                                     "to expression '{}'.".format(expression) +
                                     " Must either be numeric or 'all'.")
                if min_matching < 1:
                    raise ValueError("Invalid minimum matching parameter" +
                                     " '{}'".format(split[i]) + "passed " +
                                     "to expression '{}'.".format(expression) +
                                     " Must be > 0.")
            self.min_samples = min_matching


class FormatFilter(object):
    '''
    A class for filtering on given FORMAT fields in a VCF
    '''

    def __init__(self, vcf, expressions):
        '''
        Args:
            vcf: VariantFile object from pysam.

            expressions:
                 iterable of tuples of field names, operators and values for
                 filtering. Optional fourth parameter indicates how many
                 samples parameter has to match (either an integer or 'all' are
                 acceptable). By default only one sample passed to the 'filter'
                 method needs match an expression.

        '''
        self.vcf = vcf
        self.fields = set()
        self.expressions = []
        self._parse_expressions(expressions)

    def _parse_expressions(self, expressions):
        '''
            Args:
                expressions:
                    iterable of tuples of field names, operators and values for
                    filtering. Optional fourth parameter indicates how many
                    samples parameter has to match (either an integer or 'all'
                    are acceptable). By default only one sample passed to the
                    'filter' method needs match an expression.
        '''
        f = []
        for expression in expressions:
            self.expressions.append(FilterExpression(expression, self.vcf))

    def filter(self, record, samples):
        '''
            Read VcfRecord and returns a bitwise flag indicating
            whether each ALT allele meets parameters from self.expressions.

            record:    VcfRecord to assess according to self.expressions.

            samples:   Iterable of samples to assess.

        '''
        n_alts = len(record.alts)
        # array of bitwise flags, one flag per filter, bits set if allele meets
        # filter criteria
        flags = []
        for exp in self.expressions:
            alt_f = 0  # per-alt flag for this filter expression
            alt_counts = None
            min_smpls = 1
            if exp.min_samples < 1:
                min_smpls = len(samples)
            else:
                min_smpls = exp.min_samples
            if min_smpls > 1:
                alt_counts = [0] * n_alts
            for smp in samples:
                flt = exp.check_sample(record, smp, n_alts)
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
