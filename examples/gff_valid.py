#! /usr/local/bin/python2.7
# Copyright (C) 2014  Han Lin <hotdogee [at] gmail [dot] com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

"""
Check a GFF3 file for errors and output a validation report in markdown

Count the number of Ns in each feature, remove features with N count greater than the specified threshold. (Requires FASTA)
Check and remove features with an end coordinates larger than the landmark sequence length. (Requires FASTA or ##sequence-region)
Check if the ##sequence-region matches the FASTA file. (Requires FASTA and ##sequence-region)
Add the ##sequence-region directives if missing. (Requires FASTA)
Check and correct the phase for CDS features.

Changelog:
"""

import sys
import re
import logging
from collections import OrderedDict
from collections import defaultdict
from itertools import groupby
from urllib import quote, unquote
from textwrap import wrap
# try to import from project first
from os.path import dirname
sys.path.insert(1, dirname(dirname(__file__)))
from gff3 import Gff3

__version__ = '1.1'


def query_yes_no(question, default='yes'):
    """Ask a yes/no question via raw_input() and return their answer.

    'question' is a string that is presented to the user.
    'default' is the presumed answer if the user just hits <Enter>.
        It must be 'yes' (the default), 'no' or None (meaning
        an answer is required of the user).

    The 'answer' return value is one of 'yes' or 'no'.
    """
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('invalid default answer: "%s"' % default)

    while True:
        sys.stderr.write(question + prompt)
        choice = raw_input().strip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stderr.write('Please respond with "y" or "n".\n')
# gff_valid.py < annotations.gff > annotations.gff.validation_report
# gff_valid.py -g agla_v1_1_NALmod.gff3 > agla_v1_1_NALmod.gff3.validation_report.md
# gff_valid.py -g clec_v1_1_NALmod.gff3 > clec_v1_1_NALmod.gff3.validation_report.md
# gff_valid.py -g ofas_v1_1_NALmod.gff3 > ofas_v1_1_NALmod.gff3.validation_report.md
if __name__ == '__main__':
    logger_stderr = logging.getLogger(__name__+'stderr')
    logger_stderr.setLevel(logging.INFO)
    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger_stderr.addHandler(stderr_handler)
    logger_null = logging.getLogger(__name__+'null')
    null_handler = logging.NullHandler()
    logger_null.addHandler(null_handler)
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    Validate a GFF3 file for syntax and formating errors, parent relationship and reference sequence sanity.

    Features:
    1. Check syntax and formatting according to gff3 version 1.21.
    2. Count the number of Ns greater than the specified threshold (default: 0) in specified feature types (default: CDS). (Requires FASTA)
    3. Check for features with an end coordinates larger than the landmark sequence length. (Requires FASTA or ##sequence-region)
    4. Check if the ##sequence-region matches the FASTA file. (Requires FASTA and ##sequence-region)
    5. Check whether child features are within the coordinate boundaries of parent features.
    6. Check for the correct phase of CDS features.

    Inputs:
    1. GFF3: reads from STDIN by default, may specify the file name with the -g argument
    2. (optional) FASTA: specify the file name with the -f argument, will use the embedded ##FASTA in the GFF3 file if the external FASTA file is not specified

    Outputs:
    1. MarkDown: contains validation summary and detail sections, writes to STDOUT by default, may specify the file name with the -r argument

    Examples:
    1. Use default arguments, inout and output redirection:
        %(prog)s < a.gff > a_validation_report.txt
    2. Specify the input, output file names and options using short arguments:
        %(prog)s -g a.gff -f a.fa -n 5 -t CDS exon -r a_validation_report.txt
    3. Specify the input, output file names and options using long arguments:
        %(prog)s --gff_file a.gff --fasta_file a.fa --allowed_num_of_n 0 --check_n_feature_types CDS --report_file a_validation_report.txt
    """))
    parser.add_argument('-g', '--gff_file', type=str, help='GFF3 file to validate (default: STDIN)')
    parser.add_argument('-f', '--fasta_file', type=str, help='The external reference FASTA file for the GFF3 files, has precedence over the ##FASTA section if both exist (default: None)')
    parser.add_argument('-n', '--allowed_num_of_n', type=int, default=0,
                        help='Max number of Ns allowed in a feature, anything more will be reported as an error (default: 0)')
    parser.add_argument('-t', '--check_n_feature_types', nargs='*', default=['CDS'],
                        help='Count the number of Ns in each feature with the type specified, multiple types may be specified, ex: -t CDS exon (default: "CDS")')
    parser.add_argument('-r', '--report_file', type=str, help='Validation report file (default: STDOUT)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    test_lv = 1 # debug
    if test_lv == 0:
        args = parser.parse_args(['-g', 'annotations.gff'])
    else:
        args = parser.parse_args()

    if args.gff_file:
        logger_stderr.info('Checking GFF3 file (%s)...', args.gff_file)
    elif not sys.stdin.isatty(): # if STDIN connected to pipe or file
        args.gff_file = sys.stdin
        logger_stderr.info('Reading from STDIN...')
    else: # no input
        parser.print_help()
        sys.exit(1)

    logger_stderr.info('Checking syntax and formatting...')
    gff3 = Gff3(gff_file=args.gff_file, fasta_external=args.fasta_file, logger=logger_null)
    logger_stderr.info('Checking reference seqid, bounds and N count...')
    gff3.check_reference(allowed_num_of_n=args.allowed_num_of_n, feature_types=args.check_n_feature_types)
    logger_stderr.info('Checking parent boundaries...')
    gff3.check_parent_boundary()

    gff3.check_phase()

    if args.report_file:
        logger_stderr.info('Writing validation report (%s)...', args.report_file)
        report_fh = open(args.report_file, 'wb')
    else:
        report_fh = sys.stdout

    # Validation Summary
    report_fh.write('# GFF3 Validation Report')
    if args.gff_file and sys.stdin.isatty():
        report_fh.write(': {0:s}'.format(args.gff_file))
    report_fh.write('\n\n')

    report_fh.write('# Validation Summary\n')
    error_lines = [line for line in gff3.lines if line['line_errors']]
    if len(error_lines) == 0:
        report_fh.write('* Found 0 errors\n')
    else:
        error_list = [error for line in error_lines for error in line['line_errors']]
        error_types = sorted(list(set([error['error_type'] for error in error_list])))
        for error_type in error_types:
            report_fh.write('* Found {0:d} {1:s} errors in {2:d} lines\n'.format(
                len([error for error in error_list if error['error_type'] == error_type]), error_type,
                len([line for line in error_lines if [error for error in line['line_errors'] if error['error_type'] == error_type]])))

        report_fh.write('\n')
        report_fh.write('# Detected Errors\n')
        for line in error_lines:
            report_fh.write('* Line {0:d}: {1:s}\n'.format(line['line_index'] + 1, line['line_raw'].strip()))
            for error in line['line_errors']:
                report_fh.write('\t- {error_type}: {message}\n'.format(error_type=error['error_type'], message=error['message']))

    if args.report_file:
        report_fh.close()