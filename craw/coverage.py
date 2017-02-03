

def get_coverage(sam_file, annot_entry, start=None, stop=None, qual_thr=15, max_left=0, max_right=0):
    """
    Compute the coverage for a region position by position on each strand

    :param sam_file: the samfile openend with pysam
    :type sam_file: :class:`pysam.AlignmentFile` object.
    :param annot_entry: an entry of the annotation file
    :type annot_entry: :class:`annotation.Entry` object
    :param start: The position to start to compute the coverage(coordinates are 0-based, start position is included).
    :type start: int
    :param stop: The position to start to compute the coverage (coordinates are 0-based, stop position is excluded).
    :type stop: int
    :param qual_thr: The quality threshold
    :type qual_thr: int
    :param max_left: The highest number of base before the reference position to take in account.
    :type max_left: int
    :param max_right: The highest number of base after  the reference position to take in account.
    :type max_right: int
    :return: the coverage (all bases)
    :rtype: tuple of 2 list containing int
    """
    start < stop
    def on_forward(al_seg):
        """
        :param al_seg: a pysam aligned segment (the object used by pysam to represent an aligned read)
        :type al_seg: :class:`pysam.AlignedSegment`
        :return: True if read is mapped to forward strand
        :rtype: boolean
        """
        return not al_seg.is_reverse


    def on_reverse(al_seg):
        """
        :param al_seg: a pysam aligned segment (the object used by pysam to represent an aligned read)
        :type al_seg: :class:`pysam.AlignedSegment`
        :return: True if read is mapped to reverse strand.
        :rtype: boolean
        """
        return al_seg.is_reverse


    def coverage_one_strand(sam_file, chromosome, start, stop, qual, strand):
        """
        Compute the coverage for each position between start and stop on the chromosome on the strand.
        
        :param sam_file: the sam alignment to use
        :type sam_file: a :class:`pysam.AlignmentFile` object
        :param chromosome: the name of the chromosome
        :type chromosome: basestring
        :param start: The position to start to compute the coverage(coordinates are 0-based, start position is included).
        :type start: int
        :param stop:The position to start to compute the coverage (coordinates are 0-based, stop position is excluded).
        :type stop: int
        :param qual: The quality threshold.
        :type qual: int
        :param strand: the strand on which the read match
        :type strand: string
        :return: the coverage on forward then on reverse strand.
        The coverage is the sum of all kind bases mapped for each position
        :rtype: tuple of 2 list containing int
        """
        call_back = on_forward if strand == '+' else on_reverse
        real_start = None
        if start < 0:
            # if start is negative
            # when start is compute from large window and reads map at the beginning of the reference
            # pysam crash see issue #10
            # so we ask coverage from 0 and pad with None value for negative positions
            real_start = start
            start = 0
        try:
            coverage = sam_file.count_coverage(reference=chromosome,
                                               start=start,
                                               end=stop,
                                               quality_threshold=qual,
                                               read_callback=call_back)
        except SystemError as err:
            import sys
            print("ERROR when call count_coverage with following arguments\n",
                  "reference=", chromosome, "\n",
                  "start=", start, "\n",
                  "end=", stop, "\n",
                  "quality_threshold=", qual, "\n",
                  "read_callback=", call_back,
                  file=sys.stderr)
            raise err

        coverage = [array.tolist() for array in coverage]
        window_cov = []
        for cov_A, cov_T, cov_C, cov_G in zip(*coverage):
            window_cov.append(cov_A + cov_T + cov_C + cov_G)
        if real_start:
            window_cov = [None] * abs(real_start) + window_cov
        return window_cov

    pad_left = [None] * (max_left - (annot_entry.ref - start))
    pad_right = [None] * (max_right - (stop - annot_entry.ref))
    forward_cov = coverage_one_strand(sam_file,
                                      annot_entry.chromosome,
                                      start,
                                      stop,
                                      qual_thr,
                                      '+'
                                      )
    forward_cov = pad_left + forward_cov + pad_right
    reverse_cov = coverage_one_strand(sam_file,
                                      annot_entry.chromosome,
                                      start,
                                      stop,
                                      qual_thr,
                                      '-'
                                      )
    reverse_cov = pad_left + reverse_cov + pad_right
    return forward_cov, reverse_cov

