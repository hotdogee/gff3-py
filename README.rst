=======
gff3-py
=======

.. image:: https://github.com/hotdogee/gff3-py/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/hotdogee/gff3-py/actions/workflows/ci.yml
   :alt: CI status

.. image:: https://img.shields.io/pypi/v/gff3.svg
   :target: https://pypi.org/project/gff3/
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/pyversions/gff3.svg
   :target: https://pypi.org/project/gff3/
   :alt: Supported Python versions

.. image:: https://img.shields.io/github/license/hotdogee/gff3-py.svg
   :target: https://github.com/hotdogee/gff3-py/blob/master/LICENSE
   :alt: ISC license

Manipulate genomic features and validate the syntax and reference sequence of your |GFF3|_ files.

``gff3`` parses a |GFF3|_ file into simple Python |dict|_ and |list|_ structures,
validating the syntax as it goes and recording every problem it finds instead of
stopping at the first one. It can then check the features against their reference
sequences, from an embedded ``##FASTA`` section or an external |FASTA|_ file, and
check the phase of ``CDS`` features. You can traverse and modify the parsed features
and write the result back to a GFF3 file.

* Free software: ISC license
* Documentation: https://gff3-py.readthedocs.io
* Source code and issues: https://github.com/hotdogee/gff3-py

Features
--------

* **Simple data structures**: Parses a |GFF3|_ file into a structure composed of simple python |dict|_ and |list|_.
* **Validation**: Validates the |GFF3|_ syntax on parse, and saves the error messages in the parsed structure.
* **Best effort parsing**: Despite any detected errors, continue to parse the whole file and make as much sense to it as possible.
* Uses the python |logging|_ library to log error messages with support for custom loggers.
* Parses embedded or external |FASTA|_ sequences to check bounds and number of ``N`` s.
* Check the phase of ``CDS`` features and report the correct phase.
* Tree traversal methods ``ancestors`` and ``descendants`` return a simple ``list`` in Breadth-first search order.
* Transfer children from one feature to another using the ``adopt`` method.
* Test for overlapping features using the ``overlap`` method.
* Remove a feature and its associated features using the ``remove`` method.
* Get the sequence of a feature, or of its exons or CDSs, using the ``sequence`` method.
* Write the modified structure to a GFF3 file using the ``write`` method.

Installation
------------

``gff3`` requires Python 3.10 or newer and has no other dependencies:

.. code:: console

    $ pip install gff3

Version 1.0.1 remains available for older Python versions, including Python 2.

Quick Start
-----------

An example that just parses a GFF3 file named ``annotations.gff`` and validates it
using an external FASTA file named ``annotations.fa`` looks like:

.. code:: python

    # validate.py
    # ============
    from gff3 import Gff3

    # initialize a Gff3 object
    gff = Gff3()
    # parse GFF3 file and do syntax checking, this populates gff.lines and gff.features
    # if an embedded ##FASTA directive is found, parse the sequences into gff.fasta_embedded
    gff.parse('annotations.gff')
    # parse the external FASTA file into gff.fasta_external
    gff.parse_fasta_external('annotations.fa')
    # Check seqid, bounds and the number of Ns in each feature using one or more reference sources
    gff.check_reference(allowed_num_of_n=0, feature_types=['CDS'])
    # Checks whether child features are within the coordinate boundaries of parent features
    gff.check_parent_boundary()
    # Calculates the correct phase and checks if it matches the given phase for CDS features
    gff.check_phase()

    # every error is logged as it is found, and saved in the line it was found in
    for line in gff.lines:
        for error in line['line_errors']:
            print(line['line_index'] + 1, error['error_type'], error['message'])

Each entry of ``gff.lines`` is a ``dict`` describing one line of the file: its
``line_type`` (``'feature'``, ``'directive'``, ``'comment'`` or ``'blank'``), the raw text, the
parsed columns (``seqid``, ``type``, ``start``, ``end``, ``strand``, ``phase``,
``attributes``, ...), its ``parents`` and ``children``, and any ``line_errors``.
``gff.features`` maps each ``ID`` to the list of lines that share it. Paths can be
given as strings or ``pathlib.Path`` objects, and files given as paths are read and
written as UTF-8. If you pass an open file object instead, open it with
``encoding='utf-8', errors='surrogateescape'`` (and ``newline='\n'`` for writing)
to get the same behaviour.

The following example demonstrates how to filter, traverse, and modify the parsed gff3 ``lines`` list.

1. Change features with type ``exon`` to ``pseudogenic_exon`` and type ``transcript`` to ``pseudogenic_transcript`` if the feature has an ancestor of type ``pseudogene``

2. If a ``pseudogene`` feature overlaps with a ``gene`` feature, move all of the children from the ``pseudogene`` feature to the ``gene`` feature, and remove the ``pseudogene`` feature.

.. code:: python

    # fix_pseudogene.py
    # =================
    from gff3 import Gff3
    gff = Gff3('annotations.gff')
    type_map = {'exon': 'pseudogenic_exon', 'transcript': 'pseudogenic_transcript'}
    pseudogenes = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'pseudogene']
    for pseudogene in pseudogenes:
        # convert types
        for line in gff.descendants(pseudogene):
            if line['type'] in type_map:
                line['type'] = type_map[line['type']]
        # find overlapping gene
        overlapping_genes = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'gene' and gff.overlap(line, pseudogene)]
        if overlapping_genes:
            # move pseudogene children to overlapping gene
            gff.adopt(pseudogene, overlapping_genes[0])
            # remove pseudogene
            gff.remove(pseudogene)
    gff.write('annotations_fixed.gff')

Both scripts are in the |examples|_ directory, together with a sample ``annotations.gff``.

Validation
----------

While parsing, every line is checked against the |GFF3|_ specification, including:

* the ``##gff-version 3`` directive on the first line, and the format of the
  ``##sequence-region``, ``##species``, ``##genome-build`` and ontology directives
* nine tab separated columns, with ``.`` for empty values
* 1-based integer ``start`` and ``end`` with ``start <= end``, a numeric ``score``,
  a valid ``strand``, and a ``phase`` of 0, 1 or 2 that is required for ``CDS`` features
* percent-encoding of reserved characters in the ``seqid``, ``source``, ``type`` and attributes
* attribute syntax: one ``=`` per attribute, no repeated tags, unescaped commas,
  unknown reserved (uppercase) attributes, and the ``Target`` and ``Is_circular`` values
* duplicate ``ID`` values in non-adjacent lines and ``Parent`` references to IDs
  that are never defined

After parsing, these methods add more checks:

* ``check_reference``: the ``seqid`` exists and the feature lies within the sequence
  (using ``##sequence-region``, the embedded ``##FASTA`` or an external FASTA file),
  and the number of ``N`` s in features of the given types is not over a threshold
* ``check_parent_boundary``: child features lie within the boundaries of their parent features
* ``check_phase``: the ``phase`` of each ``CDS`` matches the phase calculated from
  the other ``CDS`` features of the same parent

Each error is a ``dict`` with an ``error_type`` (``FORMAT``, ``BOUNDS``, ``N_COUNT``,
``PHASE`` or ``STRAND``) and a ``message``.

A more feature complete GFF3 validator with a command line interface, which also generates
a validation report in Markdown, is available as |gff_valid|_:

.. code:: console

    $ python gff_valid.py -g annotations.gff -f annotations.fa -n 5 -t CDS exon -r annotations_report.md

Run ``python gff_valid.py --help`` for all of its options.

Development
-----------

Clone the repository and install it in editable mode with the development tools:

.. code:: console

    $ git clone https://github.com/hotdogee/gff3-py.git
    $ cd gff3-py
    $ python -m venv .venv
    $ source .venv/bin/activate  # on Windows: .venv\Scripts\activate
    $ python -m pip install --upgrade pip  # --group needs pip 25.1 or newer
    $ pip install -e . --group dev

Run the tests, the linter, and the build check:

.. code:: console

    $ python -m pytest
    $ ruff check
    $ python -m build
    $ twine check --strict dist/*

``tox`` runs the tests on every supported Python version. GitHub Actions runs the
tests on Linux, Windows and macOS with Python 3.10 to 3.14 for every pull request
and every push to master, and publishing a GitHub release uploads it to PyPI.

Bug reports and pull requests are welcome at https://github.com/hotdogee/gff3-py/issues.

Credits
-------

``gff3`` is written and maintained by Han Lin (`@hotdogee <https://github.com/hotdogee>`_),
with contributions from Yi Hsiao (`@hsiaoyi0504 <https://github.com/hsiaoyi0504>`_)
and Uemit Seren (`@timeu <https://github.com/timeu>`_).
See the `changelog <https://github.com/hotdogee/gff3-py/blob/master/HISTORY.rst>`_ for the release history.

.. |GFF3| replace:: ``GFF3``
.. |dict| replace:: ``dict``
.. |list| replace:: ``list``
.. |logging| replace:: ``logging``
.. |FASTA| replace:: ``FASTA``
.. |examples| replace:: ``examples/``
.. |gff_valid| replace:: ``examples/gff_valid.py``

.. _GFF3: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
.. _dict: https://docs.python.org/3/tutorial/datastructures.html#dictionaries
.. _list: https://docs.python.org/3/tutorial/datastructures.html#more-on-lists
.. _logging: https://docs.python.org/3/library/logging.html
.. _FASTA: https://en.wikipedia.org/wiki/FASTA_format
.. _examples: https://github.com/hotdogee/gff3-py/tree/master/examples
.. _gff_valid: https://github.com/hotdogee/gff3-py/blob/master/examples/gff_valid.py
