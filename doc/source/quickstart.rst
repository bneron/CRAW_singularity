.. _quickstart:

===========
Quick start
===========


craw_coverage
=============

`craw_coverage` need 2 files as inputs: the bam file and the annotation file.

 * the -b or --bam allow to specify the path to the bam file.
 * the -a --annot allow to specify the path to the annotation file.

These two options are mandatory.

.. warning::
    At the same place of `bam` file, there must be the corresponding index file (the `bam.bai` file).
    To generate the `.bai` file you have to use `samtools` program: ::

        samtools index file.bam

    see http://www.htslib.org/doc/ for more explanation.

with fix window
---------------

To compute the coverage on a fix window:
we need to specify which column name in the annotation file define the reference position.
The window will computed using this reference position.

    * --ref-col

.. note::
    if --ref-col is omitted craw_coverage will use the column position. If there not "position" column
    an error will occur.


two way to determine the window:

with --window option for a window centered on the reference position.

    * --window define the number of nucleotide to take in account before and after the reference position.

::

    craw_coverage --bam ../WTE1.bam --annot ../annotations.txt --ref-col Position --window 100

This command will compute coverage using WTE1.bam and with annotations.txt file the column used to compute the window
is 'Position' and the window length will be 100 nucleotide before the reference position and 100 nucleotides after
(201 nucleotides length).

with an asymmetrical window we have to specify two options \-\-before and \-\-after

    * \-\-before define the number of nucleotide to take in account before the reference position.
    * \-\-after  define the number of nucleotide to take in account after the reference position.

::

    craw_coverage --bam ../WTE1.bam --annot ../annotations.txt --ref-col Position --before 100 --after 500

This command will compute coverage using WTE1.bam and with annotations.txt file the column used to compute the window
is 'Position' and the window length will be 100 nucleotide before the reference position and 500 nucleotides after
(201 nucleotides length).


.. note::
    --after and --before options must be set together and are incompatible with --window option.

with variable window
--------------------

The regions must be specified in the annotation file.

* \-\-start-col define the name of the column in annotation file which define the start position of the region to compute.
* \-\-stop-col define the name of the column in annotation file which define the stop position of the region to compute.

::

    craw_coverage --bam ../WTE1.bam --annot ../annotations.txt --ref-col annotation_start --start-col annotation_start  --stop-col annotation_end

This command will compute coverage using WTE1.bam and with annotations.txt file.

* The reference position will define by the *annotation_start* column
* The first nucleotide of the window will be define by *annotation_start* column.
* The last nucleotide of the window will be define by *annotation_end* column.

other options
-------------
The folowing option are not mandatory:

* -q QUAL_THR, \-\-qual-thr QUAL_THR The minimal quality of read mapping to take it in account. (default=15)
* -s SUFFIX, \-\-suffix SUFFIX The name of the suffix to use for the output file. (default= `.cov`)
* -o OUTPUT, \-\-output OUTPUT The path of the output (default= base name of annotation file with --suffix)
* \-\-version display version information and quit.
* -h --help disply the inline help and exit.

.. warning::
    by default craw_coverage use a quality threshold of 15 (like pysam)

.. note::
    strand column mut named *strand* and can take `1/-1` or `+/-` as value for forward/reverse strands.


craw_htmp
=========


