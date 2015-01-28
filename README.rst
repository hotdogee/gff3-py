===============================
gff3-py
===============================

.. image:: https://badge.fury.io/py/gff3.png
   :target: http://badge.fury.io/py/gff3

.. image:: https://travis-ci.org/hotdogee/gff3-py.png?branch=master
   :target: https://travis-ci.org/hotdogee/gff3-py

.. image:: https://pypip.in/d/gff3/badge.png
   :target: https://pypi.python.org/pypi/gff3


Manipulate genomic features and validate the syntax and reference sequence of your |GFF3|_ files.

* Free software: BSD license
* Documentation: https://gff3-py.readthedocs.org.

Features
--------

* **Simple data structures**: Parses a |GFF3|_ file into a structure composed of simple python |dict|_ and |list|_.
* **Validation**: Validates the |GFF3|_ syntax on parse, and saves the error messages in the parsed structure.
* **Best effort parsing**: Despite any detected errors, continue to parse the whole file and make as much sense to it as possible.
* Uses the python |logging|_ library to log error messages with support for custom loggers.
* Parses embeded or external |FASTA|_ sequences to check bounds and number of ``N`` s.
* Check and correct the phase for ``CDS`` features.
* Tree traversal methods ``ancestors`` and ``descendants`` returns a simple ``list`` in Breadth-first search order.
* Transfer children and parents using the ``adopt`` and ``adopted`` methods.
* Test for overlapping features using the ``overlap`` method.
* Remove a feature and its associated features using the ``remove`` method.
* Write the modified structure to a GFF3 file using the ``write`` mthod.

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
    
A more feature complete GFF3 validator with a command line interface which also generates validation
report in MarkDown is available under ``examples/gff_valid.py``

The following example demonstrates how to filter, tranverse, and modify the parsed gff3 ``lines`` list.

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

.. |GFF3| replace:: ``GFF3``
.. |dict| replace:: ``dict``
.. |list| replace:: ``list``
.. |logging| replace:: ``logging``
.. |FASTA| replace:: ``FASTA``

.. _GFF3: http://www.sequenceontology.org/gff3.shtml
.. _dict: https://docs.python.org/2/tutorial/datastructures.html#dictionaries
.. _list: https://docs.python.org/2/tutorial/datastructures.html#more-on-lists
.. _logging: https://docs.python.org/2/library/logging.html
.. _FASTA: http://en.wikipedia.org/wiki/FASTA_format
