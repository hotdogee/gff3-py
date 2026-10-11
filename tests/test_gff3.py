#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_gff3
----------------------------------

Tests for `gff3` module.
"""

import io

import pytest

from gff3 import Gff3
from gff3.gff3 import complement, translate

from .conftest import CTG1_SEQ


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


# Sequences

def test_complement_and_translate():
    assert complement('ATGCatgcN') == 'TACGTACGN'
    assert translate('ATGGCC TAA\n') == 'MA*'
    # incomplete and unknown codons are skipped
    assert translate('ATGNNNGC') == 'M'


def test_sequence_of_a_feature(sample):
    gene = sample.features['gene1'][0]
    assert sample.sequence(gene) == CTG1_SEQ
    # a line_index works too
    assert sample.sequence(gene['line_index']) == CTG1_SEQ
    # minus strand features are reverse complemented
    match = sample.features['match1'][0]
    assert sample.sequence(match) == complement(CTG1_SEQ[1:30][::-1])
    # non-feature lines have no sequence
    assert sample.sequence(0) is None


def test_sequence_of_children(sample):
    mrna = sample.features['mRNA1'][0]
    assert sample.sequence(mrna, child_type='exon') == CTG1_SEQ[0:20] + CTG1_SEQ[30:60]
    cds = sample.sequence(mrna, child_type='CDS')
    assert cds == CTG1_SEQ[4:20] + CTG1_SEQ[30:39]
    assert translate(cds) == 'RTLASTRI'


def test_sequence_uses_the_given_reference(sample):
    reference = {'ctg1': {'id': 'ctg1', 'header': '>ctg1', 'seq': 'T' * 60}}
    assert sample.sequence(sample.features['gene1'][0], reference=reference) == 'T' * 60


def test_sequence_requires_a_reference(parse_text):
    gff = parse_text('##gff-version 3\nctg1\t.\tgene\t1\t10\t.\t+\t.\tID=g1\n')
    with pytest.raises(Exception, match='fasta reference needed'):
        gff.sequence(1)


# Parsing

HEADER = '##gff-version 3\n'


def test_parse_sample(sample):
    assert [ld['line_type'] for ld in sample.lines] == ['directive'] * 2 + ['feature'] * 8 + ['directive'] * 2
    assert [ld['line_index'] for ld in sample.lines] == list(range(12))
    assert [errors(ld) for ld in sample.lines if ld['line_errors']] == []
    gff_version, sequence_region = sample.lines[0], sample.lines[1]
    assert (gff_version['directive'], gff_version['version']) == ('##gff-version', 3)
    assert (sequence_region['directive'], sequence_region['seqid'], sequence_region['start'], sequence_region['end']) == ('##sequence-region', 'ctg1', 1, 60)
    assert [ld['directive'] for ld in sample.lines[10:]] == ['###', '##FASTA']
    assert sorted(sample.features) == ['cds1', 'exon1', 'exon2', 'gene1', 'mRNA1', 'mRNA2', 'match1']
    # lines that share an ID are one feature
    assert [ld['line_index'] for ld in sample.features['cds1']] == [7, 8]


def test_parse_feature_fields(sample):
    exon = sample.features['exon1'][0]
    assert [exon[k] for k in ('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase')] == ['ctg1', 'example', 'exon', 1, 20, '.', '+', '.']
    assert exon['attributes'] == {'ID': 'exon1', 'Parent': ['mRNA1', 'mRNA2']}
    assert exon['line_raw'] == 'ctg1\texample\texon\t1\t20\t.\t+\t.\tID=exon1;Parent=mRNA1,mRNA2\n'
    assert exon['line_status'] == 'normal'
    assert sample.features['gene1'][0]['attributes'] == {'ID': 'gene1', 'Name': 'EDEN', 'note_lowercase': 'free text'}
    assert sample.features['cds1'][1]['phase'] == 2
    match = sample.features['match1'][0]
    assert (match['score'], match['strand']) == (42.5, '-')


def test_constructor_parses_gff_and_external_fasta(sample_path, tmp_path, quiet_logger):
    fasta_path = tmp_path / 'ref.fa'
    fasta_path.write_text('>ctg1\n' + CTG1_SEQ + '\n', encoding='utf-8')
    gff = Gff3(str(sample_path), fasta_external=str(fasta_path), logger=quiet_logger)
    assert len(gff.lines) == 12
    assert gff.fasta_external['ctg1']['seq'] == CTG1_SEQ
    assert gff.fasta_embedded['ctg1']['seq'] == CTG1_SEQ


def test_parse_from_file_object(quiet_logger):
    gff = Gff3(logger=quiet_logger)
    gff.parse(io.StringIO(HEADER + 'ctg1\t.\tgene\t1\t10\t.\t+\t.\tID=g1\n'))
    assert list(gff.features) == ['g1']


def test_parse_multi_value_attributes(parse_text):
    gff = parse_text(HEADER + 'ctg1\t.\tgene\t1\t10\t.\t+\t.\tID=g1;Alias=a,b;Note=x,x;Dbxref=db:1,db:1\n')
    attributes = gff.lines[1]['attributes']
    assert attributes['Alias'] == ['a', 'b']
    # duplicate notes are allowed
    assert attributes['Note'] == ['x', 'x']
    # other duplicate values are reported and removed
    assert attributes['Dbxref'] == ['db:1']
    assert errors(gff.lines[1]) == ['Dbxref attribute has identical values (count, value): (2, db:1)']


def test_parse_target_attribute(parse_text):
    gff = parse_text(HEADER + 'ctg1\t.\tmatch\t1\t10\t.\t+\t.\tID=m1;Target=EST1 5 14 -\n')
    assert gff.lines[1]['attributes']['Target'] == {'target_id': 'EST1', 'start': 5, 'end': 14, 'strand': '-'}
    assert errors(gff.lines[1]) == []


def test_parse_directives(parse_text):
    gff = parse_text(HEADER +
                     '##feature-ontology http://example.com/so.obo\n'
                     '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7227\n'
                     '##genome-build FlyBase r6.50\n')
    assert gff.lines[1]['URI'] == 'http://example.com/so.obo'
    assert gff.lines[2]['NCBI_Taxonomy_URI'].endswith('id=7227')
    assert (gff.lines[3]['source'], gff.lines[3]['buildName']) == ('FlyBase', 'r6.50')
    assert [errors(ld) for ld in gff.lines] == [[], [], [], []]


@pytest.mark.parametrize('line, expected', [
    ('ctg1\t.\tgene\t1\t100\t.\tx\t.\tID=g1', 'Strand has illegal characters: "x"'),
    ('ctg1\t.\tCDS\t1\t100\t.\t+\t3\tID=c1', 'Phase is not 0, 1, or 2: "3"'),
    ('ctg1\t.\tCDS\t1\t100\t.\t+\tx\tID=c1', 'Phase is not a valid integer: "x"'),
    ('ctg1\t.\tCDS\t1\t100\t.\t+\t.\tID=c1', 'Phase is required for all CDS features'),
    ('ctg1\t.\tgene\t100\t1\t.\t+\t.\tID=g1', 'Start is not less than or equal to end'),
    ('ctg1\t.\tgene\t0\t100\t.\t+\t.\tID=g1', 'Start is not a valid 1-based integer coordinate: "0"'),
    ('ctg1\t.\tgene\ta\t100\t.\t+\t.\tID=g1', 'Start is not a valid integer: "a"'),
    ('ctg1\t.\tgene\t1\tb\t.\t+\t.\tID=g1', 'End is not a valid integer: "b"'),
    ('ctg1\t.\tgene\t1\t100\thigh\t+\t.\tID=g1', 'Score is not a valid floating point number: "high"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Color=red', 'Unknown reserved (uppercase) attribute: "Color"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Is_circular=yes', 'Value of Is_circular attribute is not "true": "yes"'),
    ('ctg 1\t.\tgene\t1\t100\t.\t+\t.\tID=g1', 'Seqid must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]: "ctg 1"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Name=a;Name=b', 'Found multiple attribute tags: "Name"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Name=a,b', 'Value of Name attribute contains unescaped ",": "a,b"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Note=a, b', 'Found ", " in Note attribute, possible unescaped ",": "a, b"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Name', 'Attributes must contain one and only one equal (=) sign: "Name"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;Name=', 'Empty attribute value: "Name="'),
    ('ctg1\t.\tge%ne\t1\t100\t.\t+\t.\tID=g1', 'Type must escape the percent (%) sign and any control characters: "ge%ne"'),
    ('ctg1\t\tgene\t1\t100\t.\t+\t.\tID=g1', 'Empty field: 2, must have a "."'),
    (' ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1', 'White chars not allowed at the start of a line'),
    ('ctg1\t.\tmatch\t1\t100\t.\t+\t.\tID=m1;Target=EST1 x 10 +', 'Start value of Target attribute is not a valid integer: "x"'),
    ('ctg1\t.\tgene\t1\t100\t.\t+', 'Features should contain 9 fields, got 6: '),
    ('##sequence-region ctg1 1', 'Expecting 3 fields, got 1: '),
    ('##gff-version 2', 'Version is not "3": "2"'),
    ('##not-a-directive', 'Unknown directive'),
])
def test_parse_reports_format_errors(parse_text, line, expected):
    gff = parse_text(HEADER + line + '\n')
    messages = errors(gff.lines[1])
    assert [m for m in messages if m.startswith(expected)], messages
    assert {e['error_type'] for e in gff.lines[1]['line_errors']} == {'FORMAT'}


def test_parse_reports_missing_gff_version(parse_text):
    gff = parse_text('ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n')
    assert errors(gff.lines[0]) == ['"##gff-version" missing from the first line']


def test_parse_reports_duplicate_sequence_region(parse_text):
    gff = parse_text(HEADER + '##sequence-region ctg1 1 100\n##sequence-region ctg1 1 200\n')
    assert errors(gff.lines[1]) == []
    assert errors(gff.lines[2]) == ['##sequence-region seqid: "ctg1" may only appear once']


# Reference checks

def test_check_reference_without_errors(sample):
    assert sample.check_reference() == set()
    assert [errors(ld) for ld in sample.lines if ld['line_errors']] == []


def test_check_reference_counts_ns(sample):
    assert sample.check_reference(feature_types=['gene', 'CDS']) == {2}
    error = sample.lines[2]['line_errors'][0]
    assert error['error_type'] == 'N_COUNT'
    assert error['n_segments'] == [(40, 4)]
    assert error['message'].startswith('Found 4 Ns in gene feature')


def test_check_reference_allowed_num_of_n(sample):
    assert sample.check_reference(allowed_num_of_n=4, feature_types=['gene']) == set()


@pytest.mark.parametrize('sources, expected', [
    ({}, ['End is greater than the ##sequence-region end: 60',
          'End is greater than the embedded ##FASTA sequence length: 60']),
    ({'sequence_region': True}, ['End is greater than the ##sequence-region end: 60']),
    ({'fasta_embedded': True}, ['End is greater than the embedded ##FASTA sequence length: 60']),
])
def test_check_reference_bounds(parse_text, sources, expected):
    gff = parse_text(HEADER + '##sequence-region ctg1 1 60\n'
                     'ctg1\t.\tgene\t1\t70\t.\t+\t.\tID=g1\n'
                     '##FASTA\n>ctg1\n' + CTG1_SEQ + '\n')
    assert gff.check_reference(**sources) == {2}
    assert errors(gff.lines[2]) == expected
    assert {e['error_type'] for e in gff.lines[2]['line_errors']} == {'BOUNDS'}


def test_check_reference_unknown_seqid(parse_text):
    gff = parse_text(HEADER + '##sequence-region ctg1 1 60\n'
                     'ctg2\t.\tgene\t1\t10\t.\t+\t.\tID=g1\n'
                     '##FASTA\n>ctg1\n' + CTG1_SEQ + '\n')
    assert gff.check_reference() == {2}
    assert errors(gff.lines[2]) == ['Seqid not found in any ##sequence-region: ctg2',
                                    'Seqid not found in the embedded ##FASTA: ctg2']


def test_check_reference_several_features_on_unknown_seqid(parse_text, tmp_path):
    fasta_path = tmp_path / 'ref.fa'
    fasta_path.write_text('>ctg1\n' + CTG1_SEQ + '\n', encoding='utf-8')
    gff = parse_text(HEADER + '##sequence-region ctg1 1 60\n'
                     'ctg2\t.\tgene\t1\t10\t.\t+\t.\tID=g1\n'
                     'ctg2\t.\tmRNA\t1\t10\t.\t+\t.\tID=m1;Parent=g1\n'
                     '##FASTA\n>ctg1\n' + CTG1_SEQ + '\n')
    gff.parse_fasta_external(str(fasta_path))
    # each source reports the unknown seqid once and skips the other features on it
    assert gff.check_reference() == {2}
    assert errors(gff.lines[2]) == ['Seqid not found in any ##sequence-region: ctg2',
                                    'Seqid not found in the embedded ##FASTA: ctg2',
                                    'Seqid not found in the external FASTA file: ctg2']
    assert errors(gff.lines[3]) == []


def test_check_reference_external_fasta(parse_text, tmp_path):
    fasta_path = tmp_path / 'ref.fa'
    fasta_path.write_text('>ctg1\n' + CTG1_SEQ + '\n', encoding='utf-8')
    gff = parse_text(HEADER + 'ctg1\t.\tCDS\t35\t50\t.\t+\t0\tID=c1\n')
    gff.parse_fasta_external(str(fasta_path))
    assert gff.check_reference(fasta_external=True) == {1}
    assert errors(gff.lines[1]) == ['Found 4 Ns in CDS feature of length 15 using the external FASTA, consists of 1 segment (start, length): (40, 4)']


def test_check_reference_without_sources(parse_text):
    gff = parse_text(HEADER + 'ctg1\t.\tgene\t1\t70\t.\t+\t.\tID=g1\n')
    assert gff.check_reference() == set()
    assert errors(gff.lines[1]) == []


# Tree operations

def test_descendants_and_ancestors(sample):
    gene = sample.features['gene1'][0]
    assert [ld['line_index'] for ld in sample.descendants(gene)] == [3, 4, 5, 6, 7, 8]
    assert [ld['line_index'] for ld in sample.descendants(gene['line_index'])] == [3, 4, 5, 6, 7, 8]
    assert [ld['line_index'] for ld in sample.ancestors(sample.features['exon1'][0])] == [3, 4, 2]
    assert sample.ancestors(gene) == []


def test_overlap(sample):
    gene, match, exon2 = sample.features['gene1'][0], sample.features['match1'][0], sample.features['exon2'][0]
    assert sample.overlap(gene, match)
    assert sample.overlap(match, gene)
    assert not sample.overlap(match, exon2)


def test_remove_marks_the_whole_feature_tree(sample):
    sample.remove(sample.features['exon2'][0])
    removed = [ld['line_index'] for ld in sample.lines if ld['line_status'] == 'removed']
    assert removed == [2, 3, 4, 5, 6, 7, 8]
    out = io.StringIO()
    sample.write(out, embed_fasta=False)
    assert out.getvalue() == ('##gff-version 3\n'
                              '##sequence-region ctg1 1 60\n'
                              'ctg1\texample\tmatch\t2\t30\t42.5\t-\t.\tID=match1\n'
                              '###\n')


def test_adopt_moves_children(parse_text):
    gff = parse_text(HEADER +
                     'ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=geneA\n'
                     'ctg1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=geneA\n'
                     'ctg1\t.\tgene\t1\t200\t.\t+\t.\tID=geneB\n')
    moved = gff.adopt('geneA', 'geneB')
    mrna = gff.features['mRNA1'][0]
    assert moved == [mrna]
    assert mrna['attributes']['Parent'] == ['geneB']
    assert gff.features['geneA'][0]['children'] == []
    assert gff.features['geneB'][0]['children'] == [mrna]
    assert [ld['attributes']['ID'] for ld in gff.ancestors(mrna)] == ['geneB']


def test_type_tree(sample):
    roots = sample.type_tree()
    assert [node.value for node in roots] == ['gene']
    assert repr(roots[0]) == "'gene'\n\t'mRNA'\n\t\t'CDS'\n\t\t'exon'\n"
