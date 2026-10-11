# -*- coding: utf-8 -*-
"""Run the scripts in examples/ to make sure they keep working."""
import shutil
import subprocess
import sys
from collections import Counter

import pytest

from gff3 import Gff3

from .conftest import EXAMPLES_DIR

pytestmark = pytest.mark.skipif(not (EXAMPLES_DIR / 'annotations.gff').exists(), reason='examples/ is not available')


def run_example(script, cwd, *args, **kwargs):
    return subprocess.run([sys.executable, str(EXAMPLES_DIR / script)] + list(args), cwd=str(cwd),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=300, **kwargs)


@pytest.fixture
def workdir(tmp_path):
    for name in ('annotations.gff', 'phase_test.gff3'):
        shutil.copy(str(EXAMPLES_DIR / name), str(tmp_path / name))
    return tmp_path


def test_validate(workdir):
    result = run_example('validate.py', workdir, stdin=subprocess.DEVNULL)
    assert b'PHASE: Wrong phase' in result.stderr


def test_fix_pseudogene(workdir, quiet_logger):
    run_example('fix_pseudogene.py', workdir, stdin=subprocess.DEVNULL)
    fixed = Gff3(str(workdir / 'annotations_fixed.gff'), logger=quiet_logger)
    original = Gff3(str(workdir / 'annotations.gff'), logger=quiet_logger)
    fixed_types = Counter(ld['type'] for ld in fixed.lines if ld['line_type'] == 'feature')
    original_types = Counter(ld['type'] for ld in original.lines if ld['line_type'] == 'feature')
    # the one pseudogene overlapping a gene is removed, and its transcript and exons are renamed and moved
    assert original_types - fixed_types == Counter({'pseudogene': 1, 'transcript': 1, 'exon': 4})
    assert fixed_types - original_types == Counter({'pseudogenic_transcript': 1, 'pseudogenic_exon': 4})
    transcript = [ld for ld in fixed.lines if ld['type'] == 'pseudogenic_transcript'][0]
    assert [ld['type'] for ld in fixed.ancestors(transcript)] == ['gene']


@pytest.mark.parametrize('script', ['gff_valid.py', 'gff_fix.py'])
def test_validation_report_from_file(workdir, script):
    run_example(script, workdir, '-g', 'phase_test.gff3', '-r', 'report.md', stdin=subprocess.DEVNULL)
    report = (workdir / 'report.md').read_text(encoding='utf-8')
    # the heading may include the file name, depending on whether stdin looks like a terminal
    assert report.startswith('# GFF3 Validation Report')
    assert '\n\n# Validation Summary\n* Found 13 PHASE errors in 13 lines\n\n# Detected Errors\n' in report
    assert '\t- PHASE: Wrong phase 1, should be 2\n' in report


@pytest.mark.parametrize('script', ['gff_valid.py', 'gff_fix.py'])
def test_validation_report_from_stdin(workdir, script):
    with open(str(workdir / 'annotations.gff'), 'rb') as stdin:
        result = run_example(script, workdir, stdin=stdin)
    report = result.stdout.decode('utf-8').replace('\r\n', '\n')
    assert '* Found 627 PHASE errors in 627 lines\n' in report
    assert '* Found 41 FORMAT errors in 35 lines\n' in report
