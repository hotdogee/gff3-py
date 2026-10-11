# -*- coding: utf-8 -*-
"""Tests for Gff3.write and parse -> write -> parse round trips."""
import io

import pytest

from gff3 import Gff3

from .conftest import CTG1_SEQ, EXAMPLES_DIR, SAMPLE_GFF3

FIELD_KEYS = ('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase')


def feature_records(gff):
    """Order-independent summary of every feature line that write() would output."""
    records = []
    for line_data in gff.lines:
        if line_data['line_type'] == 'feature' and line_data['line_status'] != 'removed':
            fields = tuple(line_data[k] for k in FIELD_KEYS)
            attributes = tuple(sorted((k, repr(v)) for k, v in line_data['attributes'].items()))
            records.append((fields, attributes))
    return sorted(records, key=repr)


def directive_records(gff):
    records = []
    for line_data in gff.lines:
        if line_data['line_type'] == 'directive' and line_data['directive'] not in ('###', '##FASTA'):
            if line_data['directive'] == '##sequence-region':
                records.append((line_data['directive'], line_data['seqid'], line_data['start'], line_data['end']))
            else:
                records.append((line_data['directive'], line_data['line_raw'].strip()))
    return sorted(records, key=repr)


def fasta_records(fasta_dict):
    return {key: (entry['header'], entry['seq']) for key, entry in fasta_dict.items()}


def round_trip(gff, tmp_path, quiet_logger, **write_kwargs):
    out_path = tmp_path / 'round_trip.gff3'
    gff.write(str(out_path), **write_kwargs)
    reparsed = Gff3(logger=quiet_logger)
    reparsed.parse(str(out_path))
    return reparsed, out_path


def test_parse_then_write_unmodified(sample_path, tmp_path, quiet_logger):
    """The exact scenario from issue #9: parse a path, then write to a path."""
    gff = Gff3(logger=quiet_logger)
    gff.parse(str(sample_path))
    gff.write(str(tmp_path / 'output.gff3'))
    assert (tmp_path / 'output.gff3').read_text(encoding='utf-8').startswith('##gff-version 3\n')


def test_round_trip_preserves_features(sample, tmp_path, quiet_logger):
    reparsed, _ = round_trip(sample, tmp_path, quiet_logger)
    assert feature_records(reparsed) == feature_records(sample)
    assert len(feature_records(reparsed)) == 8
    assert sorted(reparsed.features) == sorted(sample.features)
    assert directive_records(reparsed) == directive_records(sample)
    assert [ld['line_errors'] for ld in reparsed.lines if ld['line_errors']] == []


def test_round_trip_preserves_parent_child_links(sample, tmp_path, quiet_logger):
    reparsed, _ = round_trip(sample, tmp_path, quiet_logger)
    exon1 = reparsed.features['exon1'][0]
    assert sorted(f[0]['attributes']['ID'] for f in exon1['parents']) == ['mRNA1', 'mRNA2']
    gene = reparsed.features['gene1'][0]
    assert {ld['attributes']['ID'] for ld in reparsed.descendants(gene)} == {'mRNA1', 'mRNA2', 'exon1', 'exon2', 'cds1'}


def test_round_trip_embedded_fasta(sample, tmp_path, quiet_logger):
    assert fasta_records(sample.fasta_embedded) == {'ctg1': ('>ctg1 test contig', CTG1_SEQ)}
    reparsed, out_path = round_trip(sample, tmp_path, quiet_logger)
    assert fasta_records(reparsed.fasta_embedded) == fasta_records(sample.fasta_embedded)
    text = out_path.read_text(encoding='utf-8')
    assert '##FASTA\n>ctg1 test contig\n' + CTG1_SEQ + '\n' in text
    # ##sequence-region is regenerated from the FASTA sequence length
    assert '##sequence-region ctg1 1 60\n' in text


def test_round_trip_twice_is_stable(sample, tmp_path, quiet_logger):
    first, first_path = round_trip(sample, tmp_path, quiet_logger)
    second_path = tmp_path / 'second.gff3'
    first.write(str(second_path))
    assert second_path.read_bytes() == first_path.read_bytes()


def test_write_without_fasta(sample, tmp_path, quiet_logger):
    reparsed, out_path = round_trip(sample, tmp_path, quiet_logger, embed_fasta=False)
    assert '##FASTA' not in out_path.read_text(encoding='utf-8')
    assert reparsed.fasta_embedded == {}
    assert feature_records(reparsed) == feature_records(sample)


def test_write_fasta_char_limit(sample, tmp_path, quiet_logger):
    reparsed, out_path = round_trip(sample, tmp_path, quiet_logger, fasta_char_limit=25)
    fasta_part = out_path.read_text(encoding='utf-8').split('##FASTA\n')[1]
    assert fasta_part == '>ctg1 test contig\n{0}\n{1}\n{2}\n'.format(CTG1_SEQ[:25], CTG1_SEQ[25:50], CTG1_SEQ[50:])
    assert reparsed.fasta_embedded['ctg1']['seq'] == CTG1_SEQ


def test_write_to_text_file_object(sample):
    out = io.StringIO()
    sample.write(out)
    assert out.getvalue().startswith('##gff-version 3\n')
    assert '\tID=gene1;Name=EDEN;note_lowercase=free text\n' in out.getvalue()


def test_parse_and_write_accept_pathlib_paths(sample_path, tmp_path, quiet_logger):
    gff = Gff3(logger=quiet_logger)
    gff.parse(sample_path)
    out_path = tmp_path / 'pathlib.gff3'
    gff.write(out_path)
    reparsed = Gff3(out_path, logger=quiet_logger)
    assert feature_records(reparsed) == feature_records(gff)


def test_write_uses_lf_line_endings(sample, tmp_path):
    out_path = tmp_path / 'lf.gff3'
    sample.write(str(out_path))
    assert b'\r' not in out_path.read_bytes()


def test_crlf_input_round_trip(tmp_path, quiet_logger):
    crlf_path = tmp_path / 'crlf.gff3'
    crlf_path.write_bytes(SAMPLE_GFF3.replace('\n', '\r\n').encode('utf-8'))
    gff = Gff3(str(crlf_path), logger=quiet_logger)
    assert [ld['line_errors'] for ld in gff.lines if ld['line_errors']] == []
    assert gff.fasta_embedded['ctg1']['seq'] == CTG1_SEQ
    reparsed, out_path = round_trip(gff, tmp_path, quiet_logger)
    assert b'\r' not in out_path.read_bytes()
    assert feature_records(reparsed) == feature_records(gff)


def test_non_ascii_and_non_utf8_bytes_are_preserved(tmp_path, quiet_logger):
    # UTF-8 'beta' in one attribute and a stray Big5-encoded 'beta' (not valid UTF-8) in another,
    # like examples/phase_test.gff3
    g1 = b'ctg1\t.\tgene\t1\t10\t.\t+\t.\tID=g1;symbol=\xce\xb2-Spec\n'
    g2 = b'ctg1\t.\tgene\t11\t20\t.\t+\t.\tID=g2;symbol=\xa3]-Spec\n'
    in_path = tmp_path / 'bytes.gff3'
    in_path.write_bytes(b'##gff-version 3\n' + g1 + g2)
    gff = Gff3(str(in_path), logger=quiet_logger)
    assert gff.features['g1'][0]['attributes']['symbol'] == 'β-Spec'
    out_path = tmp_path / 'bytes_out.gff3'
    gff.write(str(out_path))
    assert out_path.read_bytes() == b'##gff-version 3\n' + g1 + b'###\n' + g2 + b'###\n'


def test_round_trip_example_annotations(tmp_path, quiet_logger):
    example = EXAMPLES_DIR / 'annotations.gff'
    if not example.exists():
        pytest.skip('examples/annotations.gff is not available')
    gff = Gff3(str(example), logger=quiet_logger)
    reparsed, _ = round_trip(gff, tmp_path, quiet_logger)
    assert len(feature_records(gff)) > 3000
    assert feature_records(reparsed) == feature_records(gff)
    assert directive_records(reparsed) == directive_records(gff)


@pytest.mark.parametrize('target', ['EST1 5 33 +', 'EST1 5 33 -', 'EST1 5 33'])
def test_round_trip_target_attribute(parse_text, tmp_path, quiet_logger, target):
    gff = parse_text('##gff-version 3\nctg1\t.\tmatch\t2\t30\t.\t+\t.\tID=match1;Target=' + target + '\n')
    out = io.StringIO()
    gff.write(out)
    assert '\tID=match1;Target=' + target + '\n' in out.getvalue()
    reparsed, _ = round_trip(gff, tmp_path, quiet_logger)
    assert reparsed.lines[1]['attributes']['Target'] == gff.lines[1]['attributes']['Target']
    assert reparsed.lines[1]['line_errors'] == []


def test_write_keeps_sequence_region_without_fasta(parse_text, tmp_path, quiet_logger):
    gff = parse_text('##gff-version 3\n'
                     '##sequence-region ctg1 1 1000\n'
                     '##sequence-region ctg2 1 500\n'
                     'ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n'
                     'ctg2\t.\tgene\t1\t100\t.\t+\t.\tID=g2\n')
    reparsed, out_path = round_trip(gff, tmp_path, quiet_logger)
    assert out_path.read_text(encoding='utf-8') == ('##gff-version 3\n'
                                                    '##sequence-region ctg1 1 1000\n'
                                                    'ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n'
                                                    '###\n'
                                                    '##sequence-region ctg2 1 500\n'
                                                    'ctg2\t.\tgene\t1\t100\t.\t+\t.\tID=g2\n'
                                                    '###\n')
    assert directive_records(reparsed) == directive_records(gff)


@pytest.mark.parametrize('line', [
    'ctg1\t.\tgene\t1\t100\t.\t+',
    'ctg1 . gene 1 100 . + . ID=g1',
    'ctg1\t.\tCDS\t1\t100\t.\t+\t0',
])
def test_write_keeps_malformed_feature_lines(parse_text, line):
    gff = parse_text('##gff-version 3\n' + line + '\nctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g2\n')
    assert gff.lines[1]['line_errors']
    gff.check_phase()
    out = io.StringIO()
    gff.write(out)
    assert out.getvalue() == '##gff-version 3\n' + line + '\n###\nctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g2\n###\n'


def test_write_directive_on_last_line_without_newline(parse_text):
    gff = parse_text('##gff-version 3\nctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n##species http://example.com/taxon')
    out = io.StringIO()
    gff.write(out)
    assert out.getvalue() == ('##gff-version 3\n##species http://example.com/taxon\n'
                              'ctg1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n###\n')


def test_write_sequence_region_prefers_fasta_length(sample):
    # the FASTA sequence length wins over a stale ##sequence-region directive
    sample.lines[1]['end'] = 1000
    out = io.StringIO()
    sample.write(out)
    assert '##sequence-region ctg1 1 60\n' in out.getvalue()
    assert '1000' not in out.getvalue()
