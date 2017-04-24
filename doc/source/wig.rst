.. wig:


===
wig
===

This module allow to parse wig files (wig file specifications are available here:
https://wiki.nci.nih.gov/display/tcga/wiggle+format+specification, http://genome.ucsc.edu/goldenPath/help/wiggle.html).
The wig file handle by this modules slightly differ fom de canonic specifications as
it allow to specify coverage on forward and reverse strand. If the coverage score is positive
that mean that it's on the forward strand if it's negative, it's on the reverse strand.

The WigParser and helpers
=========================

The :class:`craw.wig.WigParser` allow to parse the wig file. It read the file line by line,
test the category of the line trackLine, declarationLine or dataLine and call the right method
to parse the line and build the genome object.

The classes :class:`craw.wig.VariableChunk` and :class:`craw.wig.FixedChunk` are not keep in the final data model,
 they are just used to parsed the data lines and convert the wig file information (step, span) in coverages
 for each positions.

The data model to handle the wig information
============================================

The :class:`craw.wig.Genome` objects contains :class:`craw.wig.Chromosome` (each chromosomes ar unique and the names of chromosomes are unique).
Each chromosomes contains the coverage for the both strands.
To get the coverage for region or a position just access it with indices or slices as traditional
python list, tuple, on so on. The slicing return two lists.
The first list correspond to the coverage on this particular region for the forward strand,
the second element for the reverse strand.
By default the chromosomes are initialized with 0.0 as coverage for all positions.

All information specified in the track line are stored in the ``infos`` attribute of :class:`craw.wig.Genome` as a dict.


wig API reference
=====================

  .. automodule:: craw.wig
    :members:
    :private-members:
    :special-members:

