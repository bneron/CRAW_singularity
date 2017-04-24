.. _coverage:

========
coverage
========

*coverage* module have 3 functions

- one which is an interface to pysam.
- an other with the same api to call the wig.chromosomes
- the third is a function selector accordingly the data input, wig or bam.

get_bam_coverage
================

Get coverage from `pysam` for reference (*chromosome*) for an interval of positions, a quality on both strand.
and convert the coverage return by `pysam`. A score on each position for each base (ACGT)) in a global coverage for this
position.
It also aligned all coverages return by `pysam` on the reference position by padding left and right by ``None``.

This function is called for each entry of the annotation file.

get_wig_coverage
================

Get coverage from :class:`craw.wig.Genome` instance for reference (*chromosome*) for an interval of positions, on both strand.
The quality parameter is here just to have the same signature as get_bam_coverage but will be ignores .
It also aligned all coverages return by the Chromosome on the reference position by padding left and right by ``None``.

This function is called for each entry of the annotation file.


coverage API reference
======================

  .. automodule:: craw.coverage
    :members:
    :private-members:
    :special-members:

