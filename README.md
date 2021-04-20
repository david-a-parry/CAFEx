# CAFEx

## Case/control Allele Filtering Expressions

Filter variants in a VCF based on case/control genotype metrics.

## Installation

    pip3 install git+git://github.com/david-a-parry/CAFEx.git --user


## Usage

Detailed Usage to Follow:

```
usage: cafex [-h] [-o OUTPUT] [-t CASE [CASE ...]] [-n CONTROL [CONTROL ...]]
             [--ignore_genotypes]
             [--case_expressions EXPRESSION [EXPRESSION ...]]
             [--control_expressions EXPRESSION [EXPRESSION ...]]
             [--min_case_vaf VAF] [--max_control_vaf VAF] [--vaf_ratio RATIO]
             [--info_tag TAGNAME] [--progress_interval N] [--quiet]
             vcf

Filter case/control pair VCF based on genotype and other format fields.

positional arguments:
  vcf                   Input VCF file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output VCF/BCF file
  -t CASE [CASE ...], --case CASE [CASE ...]
                        ID of case sample(s)
  -n CONTROL [CONTROL ...], --control CONTROL [CONTROL ...]
                        ID of control sample(s)
  --ignore_genotypes    Do not use GT calls for filtering. Use this option to
                        ignore genotype calls and only filter using other
                        options (below).
  --case_expressions EXPRESSION [EXPRESSION ...]
                        Require case sample genotypes to match each of these
                        expressions. If more than one --case sample is
                        provided, default behaviour is to require only one of
                        these samples to match given expressions. For example,
                        to only retain variants with a minimum GQ of 20 in at
                        least one case sample use "GQ > 20". If you require 2
                        or more samples to match an expression add the minimum
                        number as the fourth part of your expression (e.g. "GQ
                        > 20 2"). If you require ALL cases to match an
                        expression you can use "GQ > 20 all". You may also
                        combine expressions with logical (and/or) operators.
                        Parentheses are not supported. See the README for more
                        details.
  --control_expressions EXPRESSION [EXPRESSION ...]
                        Require control sample genotypes to match each of
                        these expressions. The same rules apply as for
                        --case_expressions. For example, to only retain
                        variants with an ALT AD less than 2 in all control
                        samples use "AD < 2 all".
  --min_case_vaf VAF    Minimum VAF for case genotypes. Filter variants unless
                        at least one case sample has a VAF >= this value. VAF
                        is calculated using standard AD FORMAT fields if
                        present. Non-standard fields are supported for
                        Freebayes, Platypus, Strelka and SvABA.
  --max_control_vaf VAF
                        Maximum VAF for control genotypes. Filter variants if
                        any control sample has a VAF >= this value. VAF is
                        calculated using standard AD FORMAT fields if present.
                        Non-standard fields are supported for Freebayes,
                        Platypus, Strelka and SvABA
  --vaf_ratio RATIO     Filter variants if (case VAF) / (control VAF) is less
                        than this value. If more than one sample is provided
                        for case and/or control samples the VAF ration will be
                        calculated from the samples with the maximum VAF. VAF
                        is calcaulted as for --min_case_vaf/--max_control_vaf
                        options.
  --info_tag TAGNAME    Add bitwise flag with this tag to the INFO field of
                        your output. Bits are set in the flag indicating which
                        ALT alleles match all parameters provided by the user.
  --progress_interval N
                        Report progress every N variants. Default=100_000.
  --quiet               Suppress progress messages and only show warnings.
```
