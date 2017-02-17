.. _coverage:

========
coverage
========

*coverage* module have only one function which is an interface to pysam.
it get coverage from `pysam` for reference (*chromosome*) for an interval of positions, a quality on both strand.
and covert the coverage return by `pysam`. a score on each position for each base (ACGT)) in a global coverage for this
position.
it also aligned all coverage return by `pysam` on the reference position by padding left and right by None.

This function is called for each entry of the annotation file.


coverage API reference
======================

  .. automodule:: craw.coverage
    :members:
    :private-members:
    :special-members:

