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
