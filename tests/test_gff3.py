#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_gff3
----------------------------------

Tests for `gff3` module.
"""

import io

from gff3 import Gff3


def ids(line_datas):
    return sorted(ld['attributes']['ID'] for ld in line_datas)


def errors(line_data):
    return [error['message'] for error in line_data['line_errors']]


# Parent resolution

FORWARD_REFERENCE = """\
##gff-version 3
ctg1\t.\texon\t1\t50\t.\t+\t.\tID=exon1;Parent=mRNA1
ctg1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=gene1
ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1
###
"""


def test_forward_references_are_linked(parse_text):
    gff = parse_text(FORWARD_REFERENCE)
    exon, mrna, gene = gff.features['exon1'][0], gff.features['mRNA1'][0], gff.features['gene1'][0]
    assert exon['parents'] == [[mrna]]
    assert mrna['children'] == [exon]
    assert mrna['parents'] == [[gene]]
    assert gene['children'] == [mrna]
    assert ids(gff.descendants(gene)) == ['exon1', 'mRNA1']
    assert ids(gff.ancestors(exon)) == ['gene1', 'mRNA1']
    # forward references are allowed by GFF3, they are not errors
    assert [errors(ld) for ld in gff.lines] == [[], [], [], [], []]
    assert gff.unresolved_parents == {}


def test_forward_referenced_children_are_written(parse_text):
    gff = parse_text(FORWARD_REFERENCE)
    out = io.StringIO()
    gff.write(out)
    assert out.getvalue() == (
        '##gff-version 3\n'
        'ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n'
        'ctg1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=gene1\n'
        'ctg1\t.\texon\t1\t50\t.\t+\t.\tID=exon1;Parent=mRNA1\n'
        '###\n')


def test_unresolved_parent_is_reported(parse_text):
    gff = parse_text('##gff-version 3\n'
                     'ctg1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=missing\n')
    mrna = gff.features['mRNA1'][0]
    assert errors(mrna) == ['Parent attribute has unresolved forward reference: missing']
    assert mrna['parents'] == []
    assert 'missing' not in gff.features
    assert list(gff.unresolved_parents) == ['missing']
    # validation and writing still work
    gff.check_parent_boundary()
    out = io.StringIO()
    gff.write(out)
    assert 'ID=mRNA1;Parent=missing\n' in out.getvalue()


def test_multi_line_feature_with_interleaved_lines(parse_text):
    # exons without IDs between the lines of a discontinuous CDS
    gff = parse_text('##gff-version 3\n'
                     'ctg1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1\n'
                     'ctg1\t.\texon\t1\t30\t.\t+\t.\tParent=mRNA1\n'
                     'ctg1\t.\tCDS\t10\t30\t.\t+\t0\tID=cds1;Parent=mRNA1\n'
                     'ctg1\t.\texon\t61\t100\t.\t+\t.\tParent=mRNA1\n'
                     'ctg1\t.\tCDS\t61\t90\t.\t+\t0\tID=cds1;Parent=mRNA1\n'
                     '###\n'
                     'ctg1\t.\tgene\t200\t300\t.\t+\t.\tID=mRNA1\n')
    assert [ld['line_index'] for ld in gff.features['cds1']] == [3, 5]
    assert len(gff.features['mRNA1'][0]['children']) == 4
    assert errors(gff.lines[5]) == ['Duplicate ID: "cds1" in non-adjacent lines: 4']
    assert errors(gff.lines[7]) == ['Duplicate ID: "mRNA1" in non-adjacent lines: 2']


# Validation checks

def two_cds_mrna(strand, phases, coordinates=((1, 30), (51, 80))):
    lines = ['##gff-version 3', 'ctg1\t.\tmRNA\t1\t100\t.\t{0}\t.\tID=mRNA1'.format(strand)]
    for (start, end), phase in zip(coordinates, phases):
        lines.append('ctg1\t.\tCDS\t{0}\t{1}\t.\t{2}\t{3}\tID=cds1;Parent=mRNA1'.format(start, end, strand, phase))
    return '\n'.join(lines) + '\n'


def phase_errors(gff):
    return [(ld['line_index'], e['message']) for ld in gff.lines for e in ld['line_errors'] if e['error_type'] in ('PHASE', 'STRAND')]


def test_check_phase_accepts_correct_phases(sample):
    sample.check_phase()
    assert phase_errors(sample) == []


def test_check_phase_plus_strand(parse_text):
    # the first CDS is 31 bp long, so the second one starts with 2 bases left of a codon: phase 2
    gff = parse_text(two_cds_mrna('+', ['0', '0'], coordinates=((1, 31), (51, 80))))
    gff.check_phase()
    assert phase_errors(gff) == [(3, 'Wrong phase 0, should be 2')]


def test_check_phase_minus_strand_counts_from_the_end(parse_text):
    # on the minus strand the CDS with the largest end comes first
    gff = parse_text(two_cds_mrna('-', ['0', '0'], coordinates=((1, 30), (51, 81))))
    gff.check_phase()
    assert phase_errors(gff) == [(2, 'Wrong phase 0, should be 2')]


def test_check_phase_single_cds_must_be_phase_0(parse_text):
    gff = parse_text(two_cds_mrna('+', ['1'], coordinates=((1, 30),)))
    gff.check_phase()
    assert phase_errors(gff) == [(2, 'Wrong phase 1, should be 0')]


def test_check_phase_reports_missing_phases(parse_text):
    # issue #6: a CDS phase of "." made check_phase() raise ValueError
    gff = parse_text(two_cds_mrna('+', ['.', '.'], coordinates=((1, 31), (51, 80))))
    assert errors(gff.lines[2]) == ['Phase is required for all CDS features']
    gff.check_phase()
    assert phase_errors(gff) == [(2, 'Wrong phase ., should be 0'), (3, 'Wrong phase ., should be 2')]


def test_check_phase_inconsistent_strands(parse_text):
    text = two_cds_mrna('+', ['0', '0']).replace('+\t0\tID=cds1;Parent=mRNA1\n', '-\t0\tID=cds1;Parent=mRNA1\n', 1)
    gff = parse_text(text)
    gff.check_phase()
    assert phase_errors(gff) == [(2, 'Inconsistent CDS strand with parent: mRNA1'),
                                 (3, 'Inconsistent CDS strand with parent: mRNA1')]


def test_check_parent_boundary(parse_text):
    gff = parse_text('##gff-version 3\n'
                     'ctg1\t.\tgene\t10\t100\t.\t+\t.\tID=gene1\n'
                     'ctg1\t.\tmRNA\t10\t100\t.\t+\t.\tID=mRNA1;Parent=gene1\n'
                     'ctg1\t.\texon\t5\t50\t.\t+\t.\tID=exon1;Parent=mRNA1\n')
    gff.check_parent_boundary()
    assert errors(gff.lines[2]) == []
    assert errors(gff.lines[3]) == ['This feature is not contained within the feature boundaries of parent: mRNA1: (ctg1, 10, 100)']
    assert gff.lines[3]['line_errors'][0]['error_type'] == 'BOUNDS'


def test_checks_skip_features_without_integer_coordinates(parse_text):
    # comparing '.' with an int raised TypeError on Python 3
    gff = parse_text('##gff-version 3\n'
                     '##sequence-region ctg1 1 1000\n'
                     'ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1\n'
                     'ctg1\t.\tmRNA\t.\t100\t.\t+\t.\tID=mRNA1;Parent=gene1\n'
                     'ctg1\t.\tCDS\t.\t30\t.\t+\t0\tID=cds1;Parent=mRNA1\n'
                     'ctg1\t.\tCDS\tx\t60\t.\t+\t0\tID=cds1;Parent=mRNA1\n')
    assert gff.check_reference() == set()
    gff.check_parent_boundary()
    gff.check_phase()
    assert [errors(ld) for ld in gff.lines] == [[], [], [], [], [], ['Start is not a valid integer: "x"']]
